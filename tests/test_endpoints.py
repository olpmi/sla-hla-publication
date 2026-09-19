"""The manuscript's scientific endpoints, and the ways they are easy to misquote."""
from __future__ import annotations

import pandas as pd
import pytest

from slahla_pub import checks
from slahla_pub.paths import SOURCE_DATA


def test_all_scientific_endpoints_pass():
    failed = [c.key for c in checks.run() if not c.passed]
    assert not failed, f"failing endpoints: {failed}"


def test_subgroup_unanimity_is_over_models_not_loo_ensembles():
    """0/260 and 4/152 are model unanimity; 16.5%/52.6% is a different endpoint.

    These two get conflated because both are "how often does the subgroup agree".
    One is agreement among the five canonical models, the other agreement among
    the five leave-one-seed-out ensembles. Quoting the second as the first
    overstates subgroup stability roughly sixteenfold for class I.
    """
    pred = pd.read_csv(SOURCE_DATA / "canonical_ensemble_predictions.csv")
    for cls, expected in [("class1", 0), ("class2", 4)]:
        g = pred[pred.analysis == cls]
        assert int((g.seeds_voting_ensemble_top_subgroup == 5).sum()) == expected
        # independently, from the per-seed label string
        assert int(g.per_seed_subgroup.map(
            lambda s: len(set(str(s).split(";"))) == 1).sum()) == expected


def test_every_eplet_contrast_names_its_gap_rule():
    """Under the primary reading both class II loci show positive excess mismatch;
    under the alternative reading DRB1 reverses. A contrast without its reading
    is ambiguous, so the statistics table must carry one for every row."""
    e = pd.read_csv(SOURCE_DATA / "eplet_statistics_v6.csv")
    assert set(e.reading) == {"agree_gap_skip", "agree_gap_mismatch"}
    assert e.reading.notna().all()

    def excess(scope, reading):
        r = e[(e.analysis == "class2") & (e.scope == scope) & (e.reading == reading)].iloc[0]
        return round(-float(r.difference) * 100, 2)

    assert excess("HLA-DQB1", "agree_gap_skip") == pytest.approx(14.92, abs=0.01)
    assert excess("HLA-DRB1", "agree_gap_skip") == pytest.approx(8.99, abs=0.01)
    # the alternative reading reverses the DRB1 sign
    assert excess("HLA-DRB1", "agree_gap_mismatch") < 0


def test_permutation_p_belongs_to_the_equal_position_contrast():
    """circular_shift_p tests profile_difference, not difference."""
    e = pd.read_csv(SOURCE_DATA / "eplet_statistics_v6.csv")
    r = e[(e.analysis == "class1") & (e.scope == "ALL")
          & (e.reading == "agree_gap_skip")].iloc[0]
    assert r.p_estimator.startswith("(extreme + 1) / (shifts + 1)")
    assert float(r.profile_difference) != float(r.difference)
