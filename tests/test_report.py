"""The validation report must not overstate what a run did.

Every case here builds a synthetic verification record and renders it against a
temporary output directory, so none of it depends on a pipeline run having
happened.
"""
from __future__ import annotations

import json

import pandas as pd
import pytest

from slahla_pub import report


PANEL_COLUMNS = [
    "final_figure", "panel", "panel_name", "hla_allele", "sla_allele",
    "ce_rmsd_manifest_adopted", "ce_rmsd_recomputed",
    "rmsd_displayed_in_final_figure", "max_abs_spread_between_records_A",
    "agreement",
]


def _panel(name="class1_A", adopted=2.681, recomputed=None, displayed=2.7,
           agreement="agree"):
    return {
        "final_figure": "Figure 5", "panel": "A", "panel_name": name,
        "hla_allele": "A*02:01", "sla_allele": "SLA-1*02:01",
        "ce_rmsd_manifest_adopted": adopted, "ce_rmsd_recomputed": recomputed,
        "rmsd_displayed_in_final_figure": displayed,
        "max_abs_spread_between_records_A": 0.0001, "agreement": agreement,
    }


def _figure(name, regenerated=True, mean_diff=0.01, shift=2):
    return {"figure": name, "regenerated": regenerated, "reference": True,
            "size_match": True, "mean_abs_difference": mean_diff,
            "ink_bbox_max_shift_px": shift, "pdf_written": regenerated,
            "svg_written": regenerated}


def _record(**over):
    rec = {
        "figures": [_figure("Figure1"), _figure("Figure2")],
        "figures_missing": [], "figures_blocked": [], "blocked_outputs": [],
        "missing_inputs": [],
        "tables": {"S1": [], "S6": []},
        "recomputed_vs_archived": {"canonical_ensemble_predictions.csv": []},
        "structural_panels": {"n_displayed_panels": 1, "n_agree": 1, "disagreeing": []},
        "scientific_checks_passed": 47, "scientific_checks_total": 47,
        "source_data_unchanged": True, "source_data_changed_files": [],
        "source_data_file_count": 85,
    }
    rec.update(over)
    return rec


@pytest.fixture
def env(tmp_path, monkeypatch):
    """Point the renderer at an empty temporary output tree."""
    reports, tables = tmp_path / "reports", tmp_path / "tables"
    reports.mkdir(), tables.mkdir()
    monkeypatch.setattr(report, "REPORTS", reports)
    monkeypatch.setattr(report, "TABLES", tables)
    monkeypatch.setattr(report.runlog, "read", lambda: [])

    def write(panels=None, provenance=None):
        if panels is not None:
            pd.DataFrame(panels, columns=PANEL_COLUMNS).to_csv(
                reports / "structural_panel_audit.csv", index=False)
        if provenance is not None:
            (tables / "table_provenance.json").write_text(json.dumps(provenance))

    return write


def test_pymol_installed_but_superpositions_not_recomputed(env):
    """PyMOL being importable is not evidence that `make structures` ran.

    `make reproduce` invokes `structural --no-pymol`, so the provenance the
    pipeline recorded -- not the import -- decides what the report claims.
    """
    env(panels=[_panel(recomputed=None)],
        provenance={"S4a": "archived (run `make structures` with PyMOL to recompute)",
                    "S4b": "archived (run `make structures` with PyMOL to recompute)"})

    text = report.render(_record())

    assert "--no-pymol" in text
    assert "not** recomputed" in text or "not recomputed" in text
    assert "adopted CE RMSD" in text
    # The recomputed column is absent from the data, so it must not be tabulated,
    # and no empty measurement may be printed as a value.
    assert "recomputed CE RMSD" not in text
    assert "nan A" not in text.lower()


def test_recomputed_values_are_not_credited_to_this_run(env):
    """Values from an earlier `make structures` are shown, not claimed."""
    env(panels=[_panel(recomputed=2.638)],
        provenance={"S4a": "archived (run `make structures` with PyMOL to recompute)"})

    text = report.render(_record())

    assert "recomputed CE RMSD" in text
    assert "2.638" in text
    assert "not produced by the invocation this report covers" in text


def test_missing_required_figure_is_a_failure(env):
    """verify.main counts an unproduced, unblocked figure as a problem."""
    env(panels=[_panel()], provenance={"S4a": "archived"})
    rec = _record(figures=[_figure("Figure1"), _figure("Figure2", regenerated=False)],
                  figures_missing=["Figure2"])

    assert report.failures(rec) == ["figure Figure2: expected and not produced"]

    text = report.render(rec)
    assert "## Failures" in text
    assert "figure Figure2: expected and not produced" in text
    assert "NOT PRODUCED" in text
    # "None" must not be claimed when a figure is missing.
    assert "\n## Failures\n\nNone." not in text


def test_structural_panel_disagreement_is_a_failure(env):
    """A displayed panel that disagrees with its adopted measurement fails."""
    env(panels=[_panel(agreement="REVIEW")], provenance={"S4a": "archived"})
    rec = _record(structural_panels={"n_displayed_panels": 9, "n_agree": 8,
                                     "disagreeing": ["class1_F"]})

    assert report.failures(rec) == [
        "structural panels: 8 of 9 agree; disagreeing: class1_F"]

    text = report.render(rec)
    assert "disagreeing: class1_F" in text
    assert "REVIEW" in text
    assert "\n## Failures\n\nNone." not in text


def test_passing_record_reports_no_failures(env):
    """A clean record is the only thing that prints None."""
    env(panels=[_panel()], provenance={"S4a": "archived"})
    rec = _record()

    assert report.failures(rec) == []

    text = report.render(rec)
    assert "\n## Failures\n\nNone.\n" in text


@pytest.mark.parametrize("field, value, fragment", [
    ("tables", {"S1": ["column x differs"]}, "table S1"),
    ("recomputed_vs_archived", {"ens.csv": ["differs"]}, "intermediate ens.csv"),
    ("source_data_unchanged", False, "data/source_data was modified"),
    ("scientific_checks_passed", 45, "2 check(s) failed"),
])
def test_other_verify_criteria_reach_the_failure_summary(field, value, fragment):
    """The failure summary tracks verify.main on every criterion it counts."""
    assert any(fragment in f for f in report.failures(_record(**{field: value})))


def test_antialiasing_account_is_withheld_when_unsupported(env):
    """A large geometric shift is not explained away as glyph edges."""
    env(panels=[_panel()], provenance={"S4a": "archived"})

    small = report.render(_record(figures=[_figure("Figure1", mean_diff=0.01, shift=2)]))
    assert "consistent with text antialiasing" in small

    large = report.render(_record(figures=[_figure("Figure1", mean_diff=0.4, shift=90)]))
    assert "antialiasing" not in large
    assert "0.4" in large and "90 px" in large
