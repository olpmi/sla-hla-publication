"""HLA/SLA residue agreement at eplet-associated positions, recomputed.

Scientific starting point
-------------------------
``data/eplet_agreement/identity_eplet_agreement.parquet`` -- one row per
(SLA protein, HLA comparator, mature-protein position), carrying both readings
of a gap:

``agree_gap_skip``      gapped positions are excluded (the primary reading)
``agree_gap_mismatch``  gapped positions are counted as mismatches

Pairing is by sequence identity, not by classifier output: each SLA protein is
paired with its highest-identity eligible HLA allele. That is what makes this
analysis independent of the unstable subgroup assignment.

Three quantities that read alike in English and are not the same
---------------------------------------------------------------
``difference``          observation-weighted: every (pair, position) observation
                        counts once, so a position covered by 200 pairs counts
                        200 times. Reported with the pair bootstrap interval.
``profile_difference``  equal-position: collapse to one mean per position within
                        each locus, then average the per-locus contrasts weighted
                        by position count. Coverage cannot drive it.
``circular_shift_p``    tests ``profile_difference``, **not** ``difference``. The
                        eplet mask is rotated within each locus, preserving both
                        the autocorrelation of the profile and the clustering of
                        the eplet set.

Attaching the permutation p-value to ``difference`` would be wrong; see
``data/source_data/EPLET_STATISTIC_DEFINITIONS.md``. The Mann-Whitney column is
descriptive only: independence is violated when one HLA comparator serves up to
151 SLA proteins and hundreds of positions come from one alignment.

The eplet position mask
-----------------------
This analysis needs a set of residue positions (61 class I, 52 class II) derived
from the HLA Eplet Registry. The registry forbids redistribution of its tables
and we hold no permission covering derivatives, so the mask is **not shipped**:
supply it locally with ``python -m slahla_pub.eplet_mask --from-registry DIR``.
Without it this module produces nothing and says so. See DATA_LICENSES.md.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from . import eplet_mask, runlog
from .paths import EPLET_AGREEMENT, SOURCE_DATA, TABLES, ensure_outputs

N_BOOT = 2000
N_SHIFT = 10000
#: The seed the published statistics were generated under.
RNG_SEED = 20260913

READINGS = ("agree_gap_skip", "agree_gap_mismatch")
PRIMARY_READING = "agree_gap_skip"


def agreement() -> pd.DataFrame:
    return pd.read_parquet(EPLET_AGREEMENT / "identity_eplet_agreement.parquet")


def eplet_positions(which: str) -> set[int]:
    """The eplet-associated mature-protein positions for one class.

    Supplied locally, not shipped: see ``slahla_pub.eplet_mask`` and
    DATA_LICENSES.md for why, and for how to generate it from a registry export.
    """
    return eplet_mask.positions(which)


def _contrast(df: pd.DataFrame, col: str, mask: np.ndarray) -> float:
    return float(df.loc[mask, col].mean() - df.loc[~mask, col].mean())


def _profiles(d: pd.DataFrame, col: str, ep: set[int]):
    """Per locus: (position-mean agreement, eplet mask, weight).

    Computed once. The permutation then rolls the mask over these arrays rather
    than regrouping the frame ten thousand times.
    """
    prof = (d.groupby(["hla_locus", "position"])[col].mean()
            .rename("mean_agreement").reset_index())
    out = []
    for _, g in prof.groupby("hla_locus"):
        pos = g["position"].to_numpy(); vals = g["mean_agreement"].to_numpy()
        m = np.isin(pos, list(ep))
        if m.all() or not m.any():
            continue
        out.append((vals, m, len(pos)))
    return out


def _profile_difference(profiles, mask_fn=None) -> float:
    """Equal-position contrast, averaged over loci weighted by position count."""
    diffs, weights = [], []
    for vals, m, w in profiles:
        mm = mask_fn(m) if mask_fn else m
        diffs.append(vals[mm].mean() - vals[~mm].mean())
        weights.append(w)
    return float(np.average(diffs, weights=weights)) if diffs else float("nan")


def statistics(which: str, d: pd.DataFrame, col: str, scope: str, rng) -> dict:
    """The contrast, its pair-bootstrap interval, and the permutation test."""
    ep = eplet_positions(which)
    d = d.dropna(subset=[col]).copy()
    d[col] = d[col].astype(float)
    d["is_eplet"] = d["position"].isin(ep)
    obs = _contrast(d, col, d["is_eplet"].to_numpy())

    # 1. pair bootstrap -- the inferential unit is the (SLA, HLA) pair, because
    #    positions within a pair are dependent and HLA comparators repeat.
    pair_ids = d["sla"].to_numpy()
    uniq = np.array(sorted(set(pair_ids)))
    idx = {p: np.flatnonzero(pair_ids == p) for p in uniq}
    boots = []
    for _ in range(N_BOOT):
        pick = rng.choice(uniq, size=len(uniq), replace=True)
        sub = d.iloc[np.concatenate([idx[p] for p in pick])]
        boots.append(_contrast(sub, col, sub["is_eplet"].to_numpy()))
    lo, hi = np.percentile(boots, [2.5, 97.5])

    # 2. circular shift of the eplet mask within each locus.
    profiles = _profiles(d, col, ep)
    obs_profile = _profile_difference(profiles)
    null = []
    for _ in range(N_SHIFT):
        v = _profile_difference(profiles,
                                mask_fn=lambda m: np.roll(m, int(rng.integers(len(m)))))
        if np.isfinite(v):
            null.append(v)
    null = np.asarray(null)
    n_extreme = int((np.abs(null) >= abs(obs_profile)).sum()) if len(null) else 0
    p_uncorrected = float(n_extreme / len(null)) if len(null) else float("nan")
    # (extreme + 1) / (shifts + 1): the estimator cannot return zero, and its
    # floor 1/(N+1) is what the reported "p < ..." bound means.
    p_shift = float((n_extreme + 1) / (len(null) + 1)) if len(null) else float("nan")

    # 3. descriptive only -- independence is violated at position level.
    a = d.loc[d["is_eplet"], col]; b = d.loc[~d["is_eplet"], col]
    mw = stats.mannwhitneyu(a, b, alternative="two-sided")

    from .figure_helpers import format_permutation_p_table
    return dict(
        analysis=which, scope=scope, reading=col,
        n_sla=int(d["sla"].nunique()), n_hla=int(d["hla"].nunique()),
        n_pairs=int(d.groupby(["sla", "hla"]).ngroups),
        n_eplet_positions=int(d.loc[d["is_eplet"], "position"].nunique()),
        n_other_positions=int(d.loc[~d["is_eplet"], "position"].nunique()),
        n_eplet_observations=int(d["is_eplet"].sum()),
        n_other_observations=int((~d["is_eplet"]).sum()),
        eplet_mean=round(float(a.mean()), 4), other_mean=round(float(b.mean()), 4),
        difference=round(obs, 4),
        boot_ci_low=round(float(lo), 4), boot_ci_high=round(float(hi), 4),
        n_bootstrap=N_BOOT,
        profile_difference=round(obs_profile, 4),
        n_shifts=int(len(null)),
        n_shifts_at_least_as_extreme=n_extreme,
        circular_shift_p_uncorrected=p_uncorrected,
        circular_shift_p=p_shift,
        circular_shift_p_reported=format_permutation_p_table(p_shift, n_extreme, len(null)),
        p_estimator="(extreme + 1) / (shifts + 1); floor 1/(N+1)",
        mannwhitney_p_descriptive=float(f"{mw.pvalue:.3e}"),
        max_pairs_per_hla=int(d.groupby("hla")["sla"].nunique().max()),
        mean_pairs_per_hla=round(float(d.groupby("hla")["sla"].nunique().mean()), 2),
        primary_reading=(col == PRIMARY_READING),
    )


def figure4_values(stats_df: pd.DataFrame) -> pd.DataFrame:
    """Figure 4 shows mismatch (1 - agreement); the interval inverts with it."""
    d = stats_df[stats_df.scope == "ALL"].copy()
    d["mismatch_eplet_pct"] = (1 - d.eplet_mean) * 100
    d["mismatch_other_pct"] = (1 - d.other_mean) * 100
    d["excess_mismatch_pp"] = -d.difference * 100
    d["excess_ci_low_pp"] = -d.boot_ci_high * 100
    d["excess_ci_high_pp"] = -d.boot_ci_low * 100
    return d


def main(argv=None) -> int:
    ensure_outputs()
    try:
        eplet_mask.load()
    except eplet_mask.MissingExternalInput as exc:
        runlog.record(exc.what, exc.blocks, exc.how)
        print(f"eplet statistics: NOT RECOMPUTED - {exc.what}")
        print("  Table S3 and Figure 4 will be assembled from the archived aggregate "
              "statistics instead, and recorded as carried.")
        return 0
    agr = agreement()
    rows = []
    for which in ("class1", "class2"):
        a = agr[agr.analysis == which]
        scopes = [("ALL", a)] + [(f"HLA-{L}", a[a.hla_locus == L])
                                 for L in sorted(a.hla_locus.unique())]
        for reading in READINGS:
            for scope, sub in scopes:
                # A fresh stream per cell: the result does not depend on the order
                # the cells happen to be computed in.
                rng = np.random.default_rng(RNG_SEED)
                rows.append(statistics(which, sub, reading, scope, rng))
    d = pd.DataFrame(rows)
    order = (d.analysis + "|" + d.reading).map(
        {f"{c}|{r}": i for i, (c, r) in enumerate(
            [(c, r) for c in ("class1", "class2") for r in READINGS])})
    d = d.assign(_o=order).sort_values(["_o", "scope"],
                                       key=lambda s: s if s.name == "_o"
                                       else s.map(lambda v: "" if v == "ALL" else v)).drop(columns="_o")
    d.to_csv(TABLES / "eplet_statistics_v6.csv", index=False)
    figure4_values(d).to_csv(TABLES / "Figure4_plotted_values.csv", index=False)
    print(f"eplet statistics: {len(d)} rows "
          f"({N_BOOT} bootstrap resamples, {N_SHIFT} circular shifts, seed {RNG_SEED})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
