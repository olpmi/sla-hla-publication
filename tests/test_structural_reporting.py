"""Historical structural values must not be readable as current ones."""
from __future__ import annotations

import json

import pandas as pd
import pytest

from slahla_pub import tables
from slahla_pub.paths import ARTWORK, PANELS, SOURCE_DATA
from slahla_pub.structural import HISTORICAL_S4B_COLUMNS


def _render_report() -> dict:
    return json.loads((PANELS / "render_report.json").read_text())


def test_every_displayed_panel_prints_round1_of_its_measurement():
    """What a panel prints is round(whole-chain CE RMSD, 1). All nine agree."""
    rep = _render_report()
    man = pd.read_csv(SOURCE_DATA / "structural_panel_manifest_v6.csv")
    assert len(man) == 9
    for _, m in man.iterrows():
        ce = rep[m.panel_name]["ce_rmsd"]
        assert round(ce, 3) == pytest.approx(float(m.rmsd_adopted), abs=5e-4)
        assert round(ce, 1) == round(float(m.rmsd_adopted), 1)


def test_figure5F_is_C0702_and_prints_2_6():
    """The final legend and artwork name C*07:02; 3.4 A and C*07:01 are historical."""
    rep = _render_report()
    man = pd.read_csv(SOURCE_DATA / "structural_panel_manifest_v6.csv")
    f = man[man.panel_name == "class1_F"].iloc[0]
    assert f.hla_allele == "C*07:01" or f.hla_allele == "C*07:02"
    assert f.hla_allele == "C*07:02"
    assert round(rep["class1_F"]["ce_rmsd"], 1) == 2.6


def test_table_s4b_carries_no_submitted_version_field():
    s4b = tables.table_s4b()
    leaked = [c for c in HISTORICAL_S4B_COLUMNS if c in s4b.columns]
    assert not leaked, f"historical fields still in Table S4b: {leaked}"


def test_table_s4b_marks_the_alternative_model_as_not_a_panel():
    s4b = tables.table_s4b()
    alt = s4b[s4b.panel_name == "class1_F_alt"]
    assert len(alt) == 1
    assert alt.record_level.iloc[0] == "historical comparison (not a displayed panel)"
    assert pd.isna(alt.rmsd_displayed_in_final_figure.iloc[0])
    nine = s4b[s4b.record_level == "displayed panel"]
    assert len(nine) == 9
    assert nine.displayed_equals_round1_of_cealign.all()


def test_provenance_table_names_the_submitted_version_explicitly():
    prov = tables.table_s4b_provenance()
    assert "hla_allele_in_submitted_legend" in prov.columns
    assert "rmsd_printed_in_submitted_figure" in prov.columns
    alt = prov[prov.panel_name == "class1_F_alt"].iloc[0]
    assert "C*07:01" in alt.note and "C*07:02" in alt.note


def test_dictionary_moves_historical_fields_out_of_the_s4b_block():
    d = tables.data_dictionary()
    s4b_cols = set(d[d.table == "S4b"].column)
    assert not s4b_cols & set(HISTORICAL_S4B_COLUMNS)
    prov = d[d.table == "S4b_provenance"]
    assert len(prov) >= 6
    assert prov.description.str.contains("HISTORICAL").any()
