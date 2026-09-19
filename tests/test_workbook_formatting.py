"""Workbook styles and number formats must follow column NAMES, not positions.

Table S4b loses six columns and gains four, so positional style copying put an
integer format on `percent_identity_over_shared_positions` (75.58 rendered as 76)
and on `coverage_of_sla_chain_ATOMS` (0.646 rendered as 1). Stored values were
right throughout, so a value-only check missed it entirely; these tests assert
the *rendered* form.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd
import pytest

openpyxl = pytest.importorskip("openpyxl")

ROOT = Path(__file__).resolve().parents[1]
CORRECTED = ROOT / "outputs" / "tables" / "SLA_Supplementary_Tables_corrected.xlsx"
ORIGINAL = ROOT / "publication_artwork" / "tables" / "SLA_Supplementary_Tables.xlsx"

pytestmark = pytest.mark.skipif(not CORRECTED.exists(),
                                reason="run `make workbook` first")


def _builder():
    spec = importlib.util.spec_from_file_location(
        "bw", ROOT / "scripts" / "build_corrected_workbook.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def render(value, fmt: str) -> str:
    """What a spreadsheet shows for a stored value under a number format."""
    if value is None:
        return ""
    if isinstance(value, str):
        return value
    if fmt == "General":
        return f"{value:g}"
    if "." in fmt:
        return f"{value:.{len(fmt.split('.')[1])}f}"
    return f"{value:.0f}"


def cell(sheet: str, column: str, key: str, key_value):
    wb = openpyxl.load_workbook(CORRECTED)
    ws = wb[sheet]
    head = [c.value for c in ws[1]]
    j, k = head.index(column), head.index(key)
    for row in ws.iter_rows(min_row=2):
        if row[k].value == key_value:
            return row[j]
    raise AssertionError(f"{sheet}: no row with {key}={key_value!r}")


# ------------------------------------------------- the demonstrated failures
def test_identity_percentage_is_not_rendered_as_an_integer():
    c = cell("S4b", "percent_identity_over_shared_positions", "panel_name", "class1_F_alt")
    assert c.value == pytest.approx(75.58)
    assert render(c.value, c.number_format) == "75.5800", "75.58 must not display as 76"


def test_coverage_fraction_is_not_rendered_as_an_integer():
    c = cell("S4b", "coverage_of_sla_chain_ATOMS", "panel_name", "class1_A")
    assert 0 < c.value < 1
    assert render(c.value, c.number_format).startswith("0."), "a fraction must not display as 1"


@pytest.mark.parametrize("column", [
    "rmsd_cealign_CA", "rmsd_super_all_atom", "rmsd_align_CA", "rmsd_groove_cealign_CA",
    "coverage_of_hla_chain_CA", "coverage_of_sla_chain_CA",
    "percent_identity_over_shared_positions",
])
def test_every_s4b_measurement_keeps_decimal_places(column):
    wb = openpyxl.load_workbook(CORRECTED)
    ws = wb["S4b"]
    j = [c.value for c in ws[1]].index(column)
    for row in ws.iter_rows(min_row=2):
        c = row[j]
        if c.value is None:
            continue
        assert "." in c.number_format, f"{column} lost its decimal format"


@pytest.mark.parametrize("sheet,column", [
    ("S4b", "n_hla_CA"), ("S4b", "n_aligned_CA_pairs"), ("S4b", "n_positions_both_present"),
    ("S6", "n_pairs_covering"), ("S6", "position"), ("S4a", "n_residues"),
])
def test_counts_stay_integers(sheet, column):
    wb = openpyxl.load_workbook(CORRECTED)
    ws = wb[sheet]
    j = [c.value for c in ws[1]].index(column)
    c = ws.cell(row=2, column=j + 1)
    assert c.number_format == "0", f"{sheet}.{column} should be an integer format"


# --------------------------------------------------------- the corrections hold
def test_s6_coverage_correction_survives_the_workbook_build():
    c = cell("S6", "coverage", "position", 312)
    assert c.value == pytest.approx(1 / 259, abs=1e-12)
    assert render(c.value, c.number_format) == "0.0039"


def test_s4b_identity_correction_survives_the_workbook_build():
    c = cell("S4b", "percent_identity_over_shared_positions", "panel_name", "class1_F_alt")
    assert c.value == pytest.approx(75.58)


def test_renamed_columns_inherit_their_original_format():
    """The submitted-version fields keep the formats of the columns they replace."""
    m = _builder()
    orig = openpyxl.load_workbook(ORIGINAL)
    before = {str(c.value): f.number_format
              for c, f in zip(orig["S4b"][1], orig["S4b"][2])}
    wb = openpyxl.load_workbook(CORRECTED)
    ws = wb["S4b_provenance"]
    head = [c.value for c in ws[1]]
    for new, old in m.RENAMED.items():
        if new in head and old in before:
            j = head.index(new)
            assert ws.cell(row=2, column=j + 1).number_format == before[old], new


def test_original_workbook_is_untouched():
    orig = pd.read_excel(ORIGINAL, sheet_name="S4b")
    r = orig[orig.panel_name == "class1_F_alt"].iloc[0]
    assert float(r.percent_identity_over_shared_positions) == pytest.approx(74.92)
    assert "rmsd_reported_in_manuscript" in orig.columns
