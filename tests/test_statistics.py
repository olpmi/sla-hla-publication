"""The statistical machinery, on inputs small enough to check by hand."""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from slahla_pub import eplet_mask
from slahla_pub.eplet_stats import _contrast, _profile_difference, _profiles, eplet_positions
from slahla_pub.figure_helpers import format_permutation_p, format_permutation_p_table


def test_contrast_is_a_difference_of_means():
    d = pd.DataFrame({"v": [1.0, 1.0, 0.0, 0.0, 0.0]})
    m = np.array([True, True, False, False, False])
    assert _contrast(d, "v", m) == pytest.approx(1.0)


def test_profile_difference_weights_loci_by_position_count():
    """Equal-position contrast: a position covered by many pairs counts once."""
    d = pd.DataFrame({
        "hla_locus": ["A"] * 4 + ["B"] * 4,
        "position": [1, 2, 3, 4] * 2,
        "v": [1.0, 1.0, 0.0, 0.0, 0.8, 0.8, 0.4, 0.4],
    })
    ep = {1, 2}
    profiles = _profiles(d, "v", ep)
    assert _profile_difference(profiles) == pytest.approx((1.0 + 0.4) / 2)


def test_permutation_p_rendering_matches_both_published_formats():
    """The figure and the table annotate the same value differently; both are
    reproduced rather than harmonised."""
    assert format_permutation_p(9.999e-05, 0, 10000) == "p < 1 × 10⁻⁴"
    assert format_permutation_p(5.9994e-04, 5, 10000) == "p = 6.0 × 10⁻⁴"
    assert format_permutation_p_table(9.999e-05, 0, 10000) == "p < 1e-04"
    assert format_permutation_p_table(0.007099290, 70, 10000) == "p = 0.007099"


@pytest.mark.skipif(not eplet_mask.available(),
                    reason="the eplet position mask is not shipped; supply it locally "
                           "with `python -m slahla_pub.eplet_mask --from-registry DIR`")
def test_eplet_mask_has_the_published_position_counts():
    """When a locally supplied mask is present it must be the published set."""
    assert len(eplet_positions("class1")) == 61
    assert len(eplet_positions("class2")) == 52


def test_figure_s10_labels_come_from_the_statistics_not_from_literals():
    import inspect
    from slahla_pub import figures_final
    src = inspect.getsource(figures_final.profiles)
    assert "6.0 × 10" not in src and "1 × 10" not in src
    assert "format_permutation_p" in src
