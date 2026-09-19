"""The canonical estimator: the five-seed ensemble, recomputed from archived predictions.

Scientific starting point
-------------------------
The **archived original five-seed predictions** — ten probability tables, five
seeds for each of two classes, written by the original training runs. These are
the primary artefact underpinning every published classifier finding. Nothing
here trains, predicts, or loads a checkpoint.

The ensemble is the arithmetic mean of the predicted probability per label
across the five canonical seeds; locus mass is the sum within a locus; the top
locus is the argmax of that mass. Majority vote is computed as a *stability
descriptor*, never as the estimator.

Membership is closed
--------------------
Only ``armC_dedup_cluster`` seeds 42-46 are canonical. The archive also contains
``armA_full_random`` and ``armB_dedup_random`` (controlled split-design
comparison arms) and two legacy single-model tables. Those six are *not*
ensemble members, and
:func:`load_seed_matrix` refuses to read anything outside the canonical ten
rather than silently averaging a non-canonical arm into a published number.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from .paths import PREDICTIONS, SOURCE_DATA, TABLES, ensure_outputs

#: The canonical ensemble members: one protocol, five seeds.
CANONICAL_SEEDS: tuple[int, ...] = (42, 43, 44, 45, 46)
CLASSES = ("class1", "class2")

#: Arms present in the archive that must never enter the ensemble.
NON_CANONICAL = ("armA_full_random", "armB_dedup_random")

#: Locus display order, per class, as the published supplementary tables carry it.
#: Ordering strings below are produced by a *stable* sort over this sequence, so
#: loci that no sequence selects (argmax count 0) keep a fixed, reproducible
#: order instead of falling out of dictionary iteration.
LOCUS_DISPLAY_ORDER = {
    "class1": ("A", "B", "C"),
    "class2": ("DRB1", "DQB1", "DRB3", "DRB4", "DRB5", "DPB1"),
}

#: The SLA loci each class's published claims are computed over.
PUBLISHED_SCOPE = {
    "class1": lambda loci: ~loci.isin(["SLA-4", "SLA-5"]),
    "class2": lambda loci: loci.isin(["SLA-DRB1", "SLA-DQB1"]),
}


def run_tag(seed: int) -> str:
    return "armC_dedup_cluster" if seed == 42 else f"armC_dedup_cluster_seed{seed}"


def probs_path(which: str, tag: str):
    for bad in NON_CANONICAL:
        if bad in tag:
            raise ValueError(
                f"{tag!r} is a controlled comparison arm, not a canonical ensemble "
                f"member. The canonical ensemble is armC_dedup_cluster seeds "
                f"{CANONICAL_SEEDS}. Averaging another arm into it would change a "
                f"published number.")
    p = PREDICTIONS / f"sla_{which}__{tag}_probabilities.parquet"
    if not p.exists():
        raise FileNotFoundError(f"archived prediction table not found: {p}")
    return p


def label_space(which: str) -> tuple[list[str], np.ndarray]:
    """(labels in published order, locus of each label)."""
    labels = pd.read_csv(SOURCE_DATA / f"label_space_{which}.csv")["label"].tolist()
    loci = sorted({_top_level(l, labels) for l in labels})
    return labels, np.array([_top_level(l, labels) for l in labels])


def _loci_of(labels: list[str]) -> list[str]:
    """Recover the locus set from the label names (``B07``/``B_Other`` -> ``B``)."""
    return sorted({l.split("_")[0] if l.endswith("_Other")
                   else "".join(c for c in l if not c.isdigit()) or l
                   for l in labels})


def _top_level(label: str, labels: list[str]) -> str:
    """Collapse a subgroup label to its locus (``B07`` -> ``B``)."""
    for locus in sorted(_loci_of(labels), key=len, reverse=True):
        if label == f"{locus}_Other" or label.startswith(locus):
            return locus
    return label


def load_seed_matrix(which: str, seeds=CANONICAL_SEEDS):
    """(alleles, labels, P[seed, sequence, label], locus_of_label, sla_locus).

    Every seed's table is reindexed onto the first seed's allele order, so a row
    means the same sequence in every layer of the stack.
    """
    unexpected = set(seeds) - set(CANONICAL_SEEDS)
    if unexpected:
        raise ValueError(f"non-canonical seed(s) requested: {sorted(unexpected)}")
    labels, locus_of = label_space(which)
    base = pd.read_parquet(probs_path(which, run_tag(seeds[0])))
    alleles = base["allele"].tolist()
    loci_series = base.set_index("allele")["locus"]
    mats = []
    for s in seeds:
        d = pd.read_parquet(probs_path(which, run_tag(s))).set_index("allele")
        missing = set(alleles) - set(d.index)
        if missing:
            raise ValueError(f"{which} seed {s} is missing {len(missing)} sequence(s)")
        mats.append(d.loc[alleles, labels].to_numpy(dtype=float))
    return alleles, labels, np.stack(mats), locus_of, loci_series.loc[alleles].to_numpy()


def _locus_mass(P: np.ndarray, locus_of: np.ndarray):
    loci = sorted(set(locus_of))
    return loci, np.stack([P[..., locus_of == L].sum(axis=-1) for L in loci], axis=-1)


def ensemble(which: str, seeds=CANONICAL_SEEDS) -> pd.DataFrame:
    """Per SLA sequence: the ensemble estimate and the spread behind it."""
    alleles, labels, P, locus_of, sla_locus = load_seed_matrix(which, seeds)
    labels = np.array(labels)
    mean_p = P.mean(axis=0)
    loci, mean_mass = _locus_mass(mean_p, locus_of)
    loci = np.array(loci)
    _, seed_mass = _locus_mass(P, locus_of)

    order_lab = np.argsort(-mean_p, axis=1)
    order_loc = np.argsort(-mean_mass, axis=1)
    top_lab_i, top_loc_i = order_lab[:, 0], order_loc[:, 0]
    n = len(alleles); rows = np.arange(n)

    seed_lab = labels[P.argmax(axis=2)]
    seed_loc = loci[seed_mass.argmax(axis=2)]
    ens_lab = labels[top_lab_i]
    ens_loc = loci[top_loc_i]
    votes_loc = (seed_loc == ens_loc[None, :]).sum(axis=0)
    votes_lab = (seed_lab == ens_lab[None, :]).sum(axis=0)
    unanimous = np.array([len(set(seed_loc[:, j])) == 1 for j in range(n)])
    top_loc_by_seed = seed_mass[:, rows, top_loc_i]

    with np.errstate(divide="ignore", invalid="ignore"):
        ent = -(np.where(mean_p > 0, mean_p * np.log(mean_p), 0.0)).sum(1) / np.log(len(labels))

    out = pd.DataFrame({
        "analysis": which, "allele": alleles, "sla_locus": sla_locus,
        "ensemble_top_subgroup": ens_lab,
        "ensemble_top_subgroup_probability": mean_p[rows, top_lab_i].round(6),
        "ensemble_top_locus": ens_loc,
        "ensemble_top_locus_probability": mean_mass[rows, top_loc_i].round(6),
        "ensemble_second_locus": loci[order_loc[:, 1]],
        "ensemble_second_locus_probability": mean_mass[rows, order_loc[:, 1]].round(6),
        "ensemble_locus_margin":
            (mean_mass[rows, top_loc_i] - mean_mass[rows, order_loc[:, 1]]).round(6),
        "ensemble_normalised_entropy": ent.round(6),
        "seeds_voting_ensemble_top_locus": votes_loc,
        "seeds_voting_ensemble_top_subgroup": votes_lab,
        "n_seeds": len(seeds),
        "top_locus_probability_sd_across_seeds": top_loc_by_seed.std(axis=0, ddof=1).round(6),
        "top_locus_probability_min_across_seeds": top_loc_by_seed.min(axis=0).round(6),
        "top_locus_probability_max_across_seeds": top_loc_by_seed.max(axis=0).round(6),
        "all_seeds_agree_locus": unanimous,
        "at_least_4_of_5_agree_locus": votes_loc >= 4,
        "at_most_2_of_5_agree_locus": votes_loc <= 2,
        "per_seed_locus": [";".join(seed_loc[:, j]) for j in range(n)],
        "per_seed_subgroup": [";".join(seed_lab[:, j]) for j in range(n)],
    })
    for k, L in enumerate(loci):
        out[f"mass_{L}"] = mean_mass[:, k].round(6)
    return out


def leave_one_out(which: str, seeds=CANONICAL_SEEDS) -> pd.DataFrame:
    frames = [ensemble(which, seeds).assign(ensemble_id="full", dropped_seed="(none)")]
    for s in seeds:
        kept = tuple(x for x in seeds if x != s)
        frames.append(ensemble(which, kept).assign(ensemble_id=f"drop{s}", dropped_seed=str(s)))
    return pd.concat(frames, ignore_index=True)


def loo_sensitivity(which: str, loo: pd.DataFrame) -> dict:
    ids = [i for i in loo["ensemble_id"].unique() if i != "full"]
    wide_loc = loo.pivot(index="allele", columns="ensemble_id", values="ensemble_top_locus")
    wide_lab = loo.pivot(index="allele", columns="ensemble_id", values="ensemble_top_subgroup")
    wide_p = loo.pivot(index="allele", columns="ensemble_id",
                       values="ensemble_top_locus_probability")
    scope = PUBLISHED_SCOPE[which](loo.drop_duplicates("allele").set_index("allele")["sla_locus"])
    keep = scope[scope].index
    wl, wb, wp = wide_loc.loc[keep], wide_lab.loc[keep], wide_p.loc[keep]
    same_loc = (wl[ids].nunique(axis=1) == 1)
    same_lab = (wb[ids].nunique(axis=1) == 1)
    pairs = [(a, b) for i, a in enumerate(ids) for b in ids[i + 1:]]
    delta = (wp[ids].sub(wp["full"], axis=0)).abs()
    counts = {i: wl[i].value_counts().to_dict() for i in ["full", *ids]}
    return dict(
        analysis=which, n=len(keep),
        same_top_locus_all_loo=round(float(same_loc.mean()), 4),
        same_top_subgroup_all_loo=round(float(same_lab.mean()), 4),
        mean_pairwise_locus_agreement=round(float(np.mean(
            [(wl[a] == wl[b]).mean() for a, b in pairs])), 4),
        mean_pairwise_subgroup_agreement=round(float(np.mean(
            [(wb[a] == wb[b]).mean() for a, b in pairs])), 4),
        full_vs_loo_locus_agreement=round(float(np.mean(
            [(wl[i] == wl["full"]).mean() for i in ids])), 4),
        max_abs_change_top_locus_probability=round(float(delta.to_numpy().max()), 4),
        median_abs_change_top_locus_probability=round(float(np.median(delta.to_numpy())), 4),
        locus_counts=counts,
    )


def by_locus(loo: pd.DataFrame, which: str) -> pd.DataFrame:
    """Per-ensemble, per-SLA-locus mass and argmax counts (Figure 2C/2D)."""
    present = {c[len("mass_"):] for c in loo.columns if c.startswith("mass_")}
    order = [L for L in LOCUS_DISPLAY_ORDER[which] if L in present]
    if set(order) != present:
        raise ValueError(f"{which}: unexpected loci {sorted(present - set(order))}")
    rows = []
    for eid, g in loo.groupby("ensemble_id", sort=True):   # drop42..drop46, then full
        for scope, q in [("ALL", g)] + [(s, g[g.sla_locus == s])
                                        for s in sorted(g.sla_locus.unique())]:
            if not len(q):
                continue
            rec = {"analysis": which, "scope": scope, "ensemble_id": eid,
                   "dropped_seed": q.dropped_seed.iloc[0],
                   "n": len(q),
                   "mean_top_locus_probability": round(float(q.ensemble_top_locus_probability.mean()), 4)}
            counts = q.ensemble_top_locus.value_counts()
            means = {}
            for L in order:
                means[L] = round(float(q[f"mass_{L}"].mean()), 4)
                rec[f"mean_mass_{L}"] = means[L]
                rec[f"argmax_count_{L}"] = float(counts.get(L, 0))
            # Stable sorts over LOCUS_DISPLAY_ORDER: equal keys keep that order.
            rec["locus_order_by_mass"] = " > ".join(
                sorted(order, key=lambda L: -means[L]))
            rec["locus_order_by_argmax"] = " > ".join(
                sorted(order, key=lambda L: -counts.get(L, 0)))
            rows.append(rec)
    return pd.DataFrame(rows)


def seed_stability(which: str, seeds=CANONICAL_SEEDS) -> pd.DataFrame:
    """Per-seed argmax locus counts: the instability the ensemble exists to absorb.

    This is the *member's own* argmax over subgroup labels, collapsed to its
    locus -- the quantity each archived prediction table records in its
    ``top_label``/``top_level``/``top_probability`` columns. It is deliberately
    not the argmax of summed locus mass, which is the ensemble's rule and gives
    different counts.
    """
    rows = []
    for s in seeds:
        d = pd.read_parquet(probs_path(which, run_tag(s)))
        vc = d["top_level"].value_counts()
        rec = {"analysis": which, "seed": s, "n": len(d),
               "mean_top_probability": round(float(d["top_probability"].mean()), 4)}
        for L in ["A", "B", "C", "DPB1", "DQB1", "DRB1", "DRB3", "DRB4", "DRB5"]:
            rec[f"n_{L}"] = float(vc[L]) if L in vc else np.nan
        rows.append(rec)
    return pd.DataFrame(rows)


def validation_across_seeds() -> pd.DataFrame:
    """Held-out human accuracy per member, read from the archived run configs."""
    rows = []
    for which in CLASSES:
        for s in CANONICAL_SEEDS:
            cfg = json.loads((PREDICTIONS / f"{which}__{run_tag(s)}_run_config.json").read_text())
            rows.append({"analysis": which, "seed": s,
                         "n_train": cfg["n_train"], "n_val": cfg["n_val"],
                         "val_accuracy": round(cfg["final_val_accuracy"], 4),
                         "val_loss": round(cfg["final_val_loss"], 4),
                         "starved_labels": ",".join(cfg.get("starved_labels") or []) or np.nan,
                         "device": cfg.get("device")})
    d = pd.DataFrame(rows)
    g = d.groupby("analysis").val_accuracy
    d["val_accuracy_mean"] = g.transform("mean").round(4)
    d["val_accuracy_sd"] = g.transform(lambda x: x.std(ddof=1)).round(4)
    d["val_accuracy_min"] = g.transform("min").round(4)
    d["val_accuracy_max"] = g.transform("max").round(4)
    return d


STABILITY_COLUMNS = [
    "analysis", "allele", "sla_locus", "ensemble_top_locus",
    "seeds_voting_ensemble_top_locus", "seeds_voting_ensemble_top_subgroup",
    "top_locus_probability_sd_across_seeds", "top_locus_probability_min_across_seeds",
    "top_locus_probability_max_across_seeds", "all_seeds_agree_locus",
    "at_least_4_of_5_agree_locus", "at_most_2_of_5_agree_locus",
    "per_seed_locus", "per_seed_subgroup",
]


def main(argv=None) -> int:
    ensure_outputs()
    ens, loos, locus_rows, seedstab, sens = [], [], [], [], []
    for which in CLASSES:
        e = ensemble(which)
        loo = leave_one_out(which)
        ens.append(e)
        loos.append(loo)
        locus_rows.append(by_locus(loo, which))
        seedstab.append(seed_stability(which))
        sens.append(loo_sensitivity(which, loo))
    pred = pd.concat(ens, ignore_index=True)
    pred.to_csv(TABLES / "canonical_ensemble_predictions.csv", index=False)
    pred[STABILITY_COLUMNS].to_csv(
        TABLES / "canonical_ensemble_per_sequence_stability.csv", index=False)
    pd.concat(locus_rows, ignore_index=True).to_csv(
        TABLES / "leave_one_seed_out_by_locus.csv", index=False)
    pd.concat(seedstab, ignore_index=True).to_csv(
        TABLES / "canonical_seed_stability.csv", index=False)
    validation_across_seeds().to_csv(
        TABLES / "canonical_validation_across_seeds.csv", index=False)
    (TABLES / "leave_one_seed_out_sensitivity.json").write_text(json.dumps(sens, indent=2))
    print(f"ensemble: {len(pred)} sequences over {len(CANONICAL_SEEDS)} canonical seeds")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
