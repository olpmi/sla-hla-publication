"""Guards on the reproduction itself: closed membership, read-only inputs, no
dependency on the withdrawn analyses."""
from __future__ import annotations

import importlib

import pytest

from slahla_pub import ensemble
from slahla_pub.paths import PREDICTIONS, SOURCE_DATA


def test_only_ten_canonical_members_ship():
    got = sorted(p.name for p in PREDICTIONS.glob("*_probabilities.parquet"))
    assert len(got) == 10, got
    assert all("armC_dedup_cluster" in n for n in got)


def test_controlled_arms_are_refused_as_ensemble_members():
    """armA/armB are controlled split-design comparison arms, not ensemble
    members. Averaging one in would change a published number."""
    for bad in ensemble.NON_CANONICAL:
        with pytest.raises(ValueError, match="canonical"):
            ensemble.probs_path("class1", bad)


def test_non_canonical_seed_is_refused():
    with pytest.raises(ValueError, match="non-canonical seed"):
        ensemble.load_seed_matrix("class1", seeds=(42, 99))


@pytest.mark.parametrize("module", [
    "attribution", "ensemble_ig", "figures_attribution",
    "perturbation_pilot", "perturbation_production", "perturbation_validate",
    "run_increment", "complete_amended_cohort", "merge_increment", "refresh_v4",
    "figures_five_seed", "revision",
])
def test_withdrawn_analyses_are_absent(module):
    """The integrated-gradients and residue-replacement pipelines are withdrawn
    from the publication. Nothing here may import them, even lazily: a workflow
    that only succeeds because a cache hides a missing import is not
    self-contained."""
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module(f"slahla_pub.{module}")


def test_figure_builder_never_writes_into_its_inputs():
    import inspect
    from slahla_pub import figures_final
    src = inspect.getsource(figures_final)
    assert "to_csv(D /" not in src, "a generated CSV is being written into source_data"
    assert "D / 'panel_letter_validation.json'" not in src


def test_source_data_is_present_and_read_only_by_convention():
    assert (SOURCE_DATA / "eplet_statistics_v6.csv").exists()
    assert (SOURCE_DATA / "canonical_ensemble_predictions.csv").exists()
