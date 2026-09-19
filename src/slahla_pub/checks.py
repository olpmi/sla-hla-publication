"""The manuscript's scientific endpoints, each with its definition and source.

These are the numbers a reader would check the paper by. Every one names the
exact quantity it measures, because several of them have near-neighbours that
are easy to quote by mistake:

* Cross-seed **subgroup** unanimity (0/260 and 4/152) is agreement among the five
  canonical *models*. It is not ``same_top_subgroup_all_loo`` (0.1654 / 0.5263),
  which is agreement among the five *leave-one-seed-out ensembles* -- a different
  estimator over a different scope. Quoting the latter as the former overstates
  subgroup stability roughly sixteenfold.
* Every eplet contrast names its gap rule. Under the primary covered-residue
  reading (``agree_gap_skip``) excess eplet-associated mismatch is positive for
  both class II loci; under ``agree_gap_mismatch`` the DRB1 sign reverses. A
  contrast quoted without its reading is ambiguous.
* The permutation p-value belongs to ``profile_difference`` (the equal-position
  contrast), never to ``difference`` (the observation-weighted one).

``expected`` values are taken from the approved publication materials.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from .paths import ARTWORK, SOURCE_DATA, TABLES


@dataclass
class Check:
    key: str
    description: str
    source: str
    expected: object
    observed: object = None
    tolerance: float = 0.0
    unit: str = ""

    @property
    def passed(self) -> bool:
        if self.observed is None:
            return False
        if isinstance(self.expected, (int, float)) and not isinstance(self.expected, bool):
            return abs(float(self.observed) - float(self.expected)) <= self.tolerance
        return self.observed == self.expected


def _pred() -> pd.DataFrame:
    p = TABLES / "canonical_ensemble_predictions.csv"
    return pd.read_csv(p if p.exists() else SOURCE_DATA / "canonical_ensemble_predictions.csv")


def _stats() -> pd.DataFrame:
    p = TABLES / "eplet_statistics_v6.csv"
    return pd.read_csv(p if p.exists() else SOURCE_DATA / "eplet_statistics_v6.csv")


def _val() -> pd.DataFrame:
    p = TABLES / "canonical_validation_across_seeds.csv"
    return pd.read_csv(p if p.exists() else SOURCE_DATA / "canonical_validation_across_seeds.csv")


def run() -> list[Check]:
    out: list[Check] = []
    pred, stats, val = _pred(), _stats(), _val()

    # --- classifier: within-species performance ---------------------------------
    for cls, mean, sd in [("class1", 0.9035, 0.0168), ("class2", 0.9837, 0.0061)]:
        v = val[val.analysis == cls].val_accuracy
        out.append(Check(f"{cls}_val_accuracy_mean",
                         "mean held-out human accuracy over the five canonical models",
                         "canonical_validation_across_seeds.csv", mean,
                         round(float(v.mean()), 4), 5e-5))
        out.append(Check(f"{cls}_val_accuracy_sd",
                         "SD of held-out human accuracy across the five canonical models",
                         "canonical_validation_across_seeds.csv", sd,
                         round(float(v.std(ddof=1)), 4), 5e-5))

    # --- classifier: cross-species instability ----------------------------------
    for cls, n_sub, n_loc, n in [("class1", 0, 52, 260), ("class2", 4, 72, 152)]:
        g = pred[pred.analysis == cls]
        out.append(Check(f"{cls}_n_sequences", "SLA proteins in the classifier cohort",
                         "canonical_ensemble_predictions.csv", n, int(len(g))))
        out.append(Check(
            f"{cls}_unanimous_subgroup",
            "SLA proteins whose top HLA SUBGROUP is the same in all five canonical "
            "MODELS (not leave-one-seed-out ensembles)",
            "canonical_ensemble_predictions.csv: seeds_voting_ensemble_top_subgroup == 5",
            n_sub, int((g.seeds_voting_ensemble_top_subgroup == 5).sum())))
        out.append(Check(
            f"{cls}_unanimous_locus",
            "SLA proteins whose top HLA LOCUS is the same in all five canonical models",
            "canonical_ensemble_predictions.csv: all_seeds_agree_locus",
            n_loc, int(g.all_seeds_agree_locus.sum())))

    # The neighbouring endpoint, recorded so the two cannot be confused.
    sens_p = TABLES / "leave_one_seed_out_sensitivity.json"
    sens = json.loads((sens_p if sens_p.exists()
                       else SOURCE_DATA / "leave_one_seed_out_sensitivity.json").read_text())
    for rec, exp in zip(sens, (0.1654, 0.5263)):
        out.append(Check(
            f"{rec['analysis']}_same_top_subgroup_all_loo",
            "fraction of SLA proteins whose top subgroup agrees across the five "
            "LEAVE-ONE-SEED-OUT ENSEMBLES (a different endpoint from model unanimity)",
            "leave_one_seed_out_sensitivity.json", exp,
            rec["same_top_subgroup_all_loo"], 5e-5))

    # --- classifier: the reported estimator -------------------------------------
    c1 = pred[pred.analysis == "class1"]
    for locus, mass in [("A", 0.315), ("B", 0.485), ("C", 0.200)]:
        out.append(Check(f"class1_ensemble_mass_{locus}",
                         f"mean five-model ensemble probability mass on HLA-{locus}",
                         "canonical_ensemble_predictions.csv", mass,
                         round(float(c1[f"mass_{locus}"].mean()), 3), 5e-4))
    counts = c1.ensemble_top_locus.value_counts()
    for locus, n in [("A", 81), ("B", 160), ("C", 19)]:
        out.append(Check(f"class1_argmax_count_{locus}",
                         f"SLA proteins whose ensemble top locus is HLA-{locus} "
                         "(descriptive; carries leave-one-out ranges)",
                         "canonical_ensemble_predictions.csv", n, int(counts.get(locus, 0))))

    # --- eplet contrasts, each naming its gap rule ------------------------------
    expect = [
        ("class1", "ALL", 259, 61, 249, 0.5901, 0.7817, -0.1915, -0.1968, -0.1869, -0.1783),
        ("class2", "ALL", 151, 52, 165, 0.7233, 0.8367, -0.1134, -0.1223, -0.1043, -0.1366),
        ("class2", "HLA-DQB1", 53, 52, 165, 0.6979, 0.8470, -0.1492, -0.1599, -0.1377, None),
        ("class2", "HLA-DRB1", 98, 52, 165, 0.7392, 0.8291, -0.0899, -0.1009, -0.0782, None),
    ]
    for cls, scope, npairs, nep, noth, em, om, diff, lo, hi, prof in expect:
        r = stats[(stats.analysis == cls) & (stats.scope == scope)
                  & (stats.reading == "agree_gap_skip")].iloc[0]
        tag = f"{cls}_{scope.replace('-', '')}_gap_skip"
        out.append(Check(f"{tag}_n_pairs", f"{cls} {scope} identity-selected pairs "
                         "(primary covered-residue reading)",
                         "eplet_statistics_v6.csv", npairs, int(r.n_pairs)))
        out.append(Check(f"{tag}_n_eplet_positions", "eplet-associated positions",
                         "eplet_statistics_v6.csv", nep, int(r.n_eplet_positions)))
        out.append(Check(f"{tag}_n_other_positions", "other aligned positions",
                         "eplet_statistics_v6.csv", noth, int(r.n_other_positions)))
        out.append(Check(f"{tag}_difference",
                         "observation-weighted agreement contrast, gaps EXCLUDED",
                         "eplet_statistics_v6.csv", diff, float(r.difference), 5e-5))
        out.append(Check(f"{tag}_excess_mismatch_pp",
                         "excess eplet-associated MISMATCH, gaps EXCLUDED",
                         "eplet_statistics_v6.csv", round(-diff * 100, 2),
                         round(-float(r.difference) * 100, 2), 5e-3, "percentage points"))
        out.append(Check(f"{tag}_boot_ci", "pair-bootstrap 95% interval on the contrast",
                         "eplet_statistics_v6.csv", [lo, hi],
                         [float(r.boot_ci_low), float(r.boot_ci_high)]))
        if prof is not None:
            out.append(Check(f"{tag}_profile_difference",
                             "equal-position contrast; THIS is what circular_shift_p tests",
                             "eplet_statistics_v6.csv", prof, float(r.profile_difference), 5e-5))

    # The alternative reading, where the class II DRB1 sign reverses.
    r = stats[(stats.analysis == "class2") & (stats.scope == "HLA-DRB1")
              & (stats.reading == "agree_gap_mismatch")].iloc[0]
    out.append(Check("class2_HLADRB1_gap_mismatch_difference",
                     "class II DRB1 contrast with gaps COUNTED AS MISMATCHES -- "
                     "positive, i.e. the sign reverses against the primary reading",
                     "eplet_statistics_v6.csv", 0.1007, float(r.difference), 5e-5))

    # --- structural -------------------------------------------------------------
    s4a = pd.read_csv(ARTWORK / "tables" / "TableS4a_structural_model_confidence_per_model.csv")
    s4b = pd.read_csv(ARTWORK / "tables" / "TableS4b_structural_superpositions_per_pair.csv")
    out.append(Check("structural_models_scored", "individual models in Table S4a",
                     "TableS4a", 23, int(len(s4a))))
    out.append(Check("structural_pairs_scored", "model pairs in Table S4b "
                     "(nine published panels plus the Figure 5F allele alternative)",
                     "TableS4b", 10, int(len(s4b))))
    return out


def main(argv=None) -> int:
    checks = run()
    width = max(len(c.key) for c in checks)
    failed = 0
    for c in checks:
        ok = c.passed
        failed += not ok
        print(f"{'PASS' if ok else 'FAIL'}  {c.key:<{width}}  expected {c.expected!r:>28}  "
              f"observed {c.observed!r}")
    print(f"\n{len(checks) - failed}/{len(checks)} scientific endpoint checks passed")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
