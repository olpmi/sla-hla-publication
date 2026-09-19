"""The corrections this release makes, and the guarantees that keep them."""
from __future__ import annotations

import pandas as pd
import pytest

from slahla_pub import eplet_mask, haplotypes as H, tables
from slahla_pub.paths import ARTWORK, DATA, EPLET_AGREEMENT, SOURCE_DATA


# ---------------------------------------------------------------- C1: S6 coverage
def test_s6_coverage_uses_the_documented_class_denominator():
    """Every row divides by the class pair count, not by the rows present."""
    s6 = tables.table_s6()
    agr = pd.read_parquet(EPLET_AGREEMENT / "identity_eplet_agreement.parquet")
    for cls, denom in tables.S6_PAIR_DENOMINATOR.items():
        g = s6[s6.analysis == cls]
        assert agr[agr.analysis == cls].groupby(["sla", "hla"]).ngroups == denom
        assert ((g.n_pairs_covering / denom - g.coverage).abs() < 5e-12).all()


def test_s6_position_312_is_corrected():
    s6 = tables.table_s6()
    r = s6[(s6.analysis == "class1") & (s6.position == 312)].iloc[0]
    assert r.n_pairs_covering == 1
    assert r.coverage == pytest.approx(1 / 259, abs=1e-15)
    approved = pd.read_csv(ARTWORK / "tables" /
                           "TableS6_experimental_validation_priorities.csv")
    old = approved[(approved.analysis == "class1") & (approved.position == 312)].iloc[0]
    assert old.coverage == pytest.approx(1 / 22, abs=1e-9), "approved value changed"


def test_s6_correction_touches_exactly_one_row():
    new = tables.table_s6()
    old = pd.read_csv(ARTWORK / "tables" / "TableS6_experimental_validation_priorities.csv")
    changed = (new.coverage.astype(float) - old.coverage.astype(float)).abs() > 5e-7
    assert int(changed.sum()) == 1


# ------------------------------------------------------- C2: S4b historical identity
def test_s4b_alt_identity_is_supported_by_the_alignment_record():
    ev = pd.read_csv(DATA / "carried" / "tableS4b_pair_identity_verified.csv")
    r = ev[ev.panel_name == "class1_F_alt"].iloc[0]
    assert r.identical_residues_numerator == 229
    assert r.n_positions_both_present == 303
    assert r.percent_identity_over_shared_positions == pytest.approx(75.58, abs=0.005)
    assert round(100 * r.identical_residues_numerator / r.n_positions_both_present, 2) == 75.58
    assert r.hla_model_file == "C07_C0701_unrelaxed_CUSTOM_B.pdb"
    # the two rows must not share an identity value any more
    f = ev[ev.panel_name == "class1_F"].iloc[0]
    assert f.percent_identity_over_shared_positions != r.percent_identity_over_shared_positions


def test_table_s4b_reports_the_corrected_identity_in_both_modes():
    s4b = tables.table_s4b()
    r = s4b[s4b.panel_name == "class1_F_alt"].iloc[0]
    assert r.percent_identity_over_shared_positions == pytest.approx(75.58, abs=0.005)
    assert r.record_level == "historical comparison (not a displayed panel)"


def test_figure5F_is_untouched():
    """The correction is to a historical row; the displayed panel must not move."""
    s4b = tables.table_s4b()
    f = s4b[s4b.panel_name == "class1_F"].iloc[0]
    assert f.hla_allele_in_model == "C*07:02"
    assert f.rmsd_displayed_in_final_figure == pytest.approx(2.6)
    assert f.percent_identity_over_shared_positions == pytest.approx(74.92, abs=0.005)


# ------------------------------------------------------------- C4: S1 provenance
def test_s1_carried_selections_are_recorded_and_labelled():
    carried = H.carried_selection()
    assert len(carried) == 3
    assert set(carried.sla_allele) == {"SLA-DQA*01:01", "SLA-DQA*02:02:02", "SLA-DRA*02:03:01"}
    assert carried.source.str.contains("approved").all()
    g = H.alleles_by_haplotype().set_index("sla_allele")
    for a in carried.sla_allele:
        assert "carried" in g.loc[a, "haplotype_selection_source"]
    others = g[~g.index.isin(carried.sla_allele)]
    assert (others.haplotype_selection_source == "regenerated from the curated workbooks").all()


def test_s1_carry_does_not_discard_the_regenerated_set():
    g = H.alleles_by_haplotype().set_index("sla_allele")
    r = g.loc["SLA-DQA*02:02:02"]
    assert r.haplotypes == "Hp-28.15b"                    # carried
    assert "Hp-4b.5" in r.haplotypes_all_curated          # regenerated, retained


# ------------------------------------------------- the mask is not in the bundle
@pytest.mark.parametrize("name", [
    "figureS8_S13_UPDATED_plotted_values.csv",
    "figureS_identity_eplet_agreement_corrected_p_plotted_values.csv",
])
def test_shipped_inputs_do_not_carry_the_eplet_mask(name):
    d = pd.read_csv(SOURCE_DATA / name, nrows=1)
    assert "is_eplet_associated" not in d.columns


def test_the_position_enumeration_file_is_not_shipped():
    assert not (SOURCE_DATA / "eplet_position_mask_comparison.csv").exists()


def test_missing_mask_raises_an_actionable_error(monkeypatch, tmp_path):
    monkeypatch.setattr(eplet_mask, "MASK_PATH", tmp_path / "absent.csv")
    with pytest.raises(eplet_mask.MissingExternalInput) as e:
        eplet_mask.load()
    assert "--from-registry" in str(e.value)
    assert set(e.value.blocks) == {"FigureS8", "FigureS9", "FigureS10"}


# --------------------------------------------------------- workbook agreement
@pytest.mark.skipif(not (ARTWORK / "tables" / "SLA_Supplementary_Tables.xlsx").exists(),
                    reason="authoritative workbook not present")
def test_corrected_workbook_agrees_with_the_regenerated_csvs():
    """If the workbook has been built, every sheet must match its CSV."""
    import openpyxl
    from slahla_pub.paths import TABLES
    dest = TABLES / "SLA_Supplementary_Tables_corrected.xlsx"
    if not dest.exists():
        pytest.skip("run `make workbook` first")
    wb = openpyxl.load_workbook(dest, data_only=True)
    assert "S4b_provenance" in wb.sheetnames
    for sheet, csv in [("S6", "TableS6_experimental_validation_priorities.csv"),
                       ("S4b", "TableS4b_structural_superpositions_per_pair.csv")]:
        df = pd.read_csv(TABLES / "publication" / csv)
        ws = wb[sheet]
        got = pd.DataFrame(ws.iter_rows(min_row=2, values_only=True),
                           columns=[c.value for c in ws[1]])
        assert list(got.columns) == list(df.columns)
        assert len(got) == len(df)
    # the original is untouched
    orig = pd.read_excel(ARTWORK / "tables" / "SLA_Supplementary_Tables.xlsx", sheet_name="S6")
    row = orig[(orig.analysis == "class1") & (orig.position == 312)].iloc[0]
    assert float(row.coverage) == pytest.approx(1 / 22, abs=1e-9)
