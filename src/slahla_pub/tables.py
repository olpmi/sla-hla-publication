"""Assemble the publication tables S1-S6 (S4 as separate S4a and S4b).

Each table declares where it starts from. Three cases, and the code says which
applies rather than presenting a copy as a regeneration:

``recomputed``    every column is recalculated here from archived scientific
                  inputs (S2, S3, S4a, S4b, S5).
``partial``       some columns are recalculated and the rest require an upstream
                  stage this repository does not ship inputs for (S1).
``archived``      carried from the approved publication materials, because its
                  inputs cannot be redistributed (S6).

``PUBLICATION_OUTPUTS.csv`` carries the same classification per output, and
``outputs/reports/validation_report.md`` records what was compared.
"""
from __future__ import annotations

import json

import pandas as pd

from . import haplotypes as H
from .paths import (ARTWORK, DATA, EPLET_AGREEMENT, SOURCE_DATA, TABLES,
                    ensure_outputs)

NAMES = {
    "S1": "TableS1_haplotype_nearest_and_most_distant_HLA.csv",
    "S2": "TableS2_canonical_split_leakage_and_five_seed_validation.csv",
    "S3": "TableS3_HLA_eplet_associated_position_statistics.csv",
    "S4a": "TableS4a_structural_model_confidence_per_model.csv",
    "S4b": "TableS4b_structural_superpositions_per_pair.csv",
    "S4b_provenance": "TableS4b_provenance_submitted_version.csv",
    "S5": "TableS5_structural_position_annotations.csv",
    "S6": "TableS6_experimental_validation_priorities.csv",
}


def _published(key: str) -> pd.DataFrame:
    return pd.read_csv(ARTWORK / "tables" / NAMES[key])


def table_s1() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Nearest and most distant CIWD 3.0 HLA allele per SLA allele.

    Partial. The haplotype half -- which SLA allele sits in which haplotype, herd
    and local line -- is regenerated here from the curated workbooks in
    ``data/haplotypes/``. The ``n_ciwd_hits``, ``nearest_*`` and
    ``most_distant_*`` columns are pairwise sequence identities against the CIWD
    3.0 common-and-well-documented allele set, which needs the DIAMOND all-vs-all
    stage over IPD-IMGT/HLA and the CIWD release. Neither can be redistributed
    here, so those columns are carried from the approved table and the upstream
    route is documented in REPRODUCIBILITY.md.
    """
    published = _published("S1")
    recomputed = H.alleles_by_haplotype()
    return published, recomputed


def table_s2() -> pd.DataFrame:
    """Canonical split design, leakage audit and five-seed held-out validation."""
    split = pd.read_csv(SOURCE_DATA / "canonical_split_summary.csv")
    val = pd.read_csv(TABLES / "canonical_validation_across_seeds.csv")
    canon = split[split.arm == "armC_dedup_cluster"].copy()
    keep = ["analysis", "arm", "description", "dedup", "split_strategy", "n_records",
            "n_distinct_sequences", "n_clusters_99", "n_train", "n_validation",
            "actual_val_fraction", "n_labels", "labels_in_validation",
            "labels_absent_from_validation", "exact_sequence_val_share",
            "cluster_99_val_share", "two_field_val_share"]
    canon = canon[keep].rename(columns={"n_train": "n_train_split"})
    out = val.merge(canon, on="analysis", how="left").rename(
        columns={"val_accuracy_sd": "val_accuracy_std"})
    order = ["analysis", "seed", "n_train", "n_val", "val_accuracy", "val_loss",
             "starved_labels", "device", "arm", "description", "dedup", "split_strategy",
             "n_records", "n_distinct_sequences", "n_clusters_99", "n_train_split",
             "n_validation", "actual_val_fraction", "n_labels", "labels_in_validation",
             "labels_absent_from_validation", "exact_sequence_val_share",
             "cluster_99_val_share", "two_field_val_share", "val_accuracy_mean",
             "val_accuracy_std", "val_accuracy_min", "val_accuracy_max"]
    return out[[c for c in order if c in out.columns]]


def table_s3() -> pd.DataFrame:
    """HLA eplet-associated position statistics.

    Recomputed by ``eplet_stats`` when the eplet position mask is available
    locally; otherwise carried from the archived statistics, which hold group
    means, intervals and counts but no position information.
    """
    p = TABLES / "eplet_statistics_v6.csv"
    d = pd.read_csv(p if p.exists() else SOURCE_DATA / "eplet_statistics_v6.csv")
    return d.drop(columns=[c for c in ["circular_shift_p_uncorrected"] if c in d.columns])


def _use_recomputed(path, required: list[str]) -> pd.DataFrame | None:
    """The recomputed frame, but only if it is complete.

    The structural recompute needs PyMOL. Without it the default workflow carries
    the approved table instead of a partially-filled one, and says so, rather
    than silently dropping the columns PyMOL would have produced.
    """
    if not path.exists():
        return None
    d = pd.read_csv(path)
    missing = [c for c in required if c not in d.columns]
    return None if missing else d


def table_s4a() -> pd.DataFrame:
    """Per-model confidence and superposition against the experimental reference.

    Recomputed when the optional structural workflow has run (``make structures``),
    otherwise carried from the approved table.
    """
    pubtab = _published("S4a")
    cols = list(pubtab.columns)
    d = _use_recomputed(TABLES / "structural_model_confidence_recomputed.csv", cols)
    if d is None:
        return pubtab
    d = d[d.allele.notna()]
    # Published row order: HLA before SLA, class I before class II, then allele.
    d = d.sort_values(["system", "mhc_class", "allele"]).reset_index(drop=True)
    return d[[c for c in cols if c in d.columns]]


#: Alignment evidence for Table S4b's two identity columns: model checksums,
#: aligner version, numerator and denominator for every pair.
VERIFIED_IDENTITY = DATA / "carried" / "tableS4b_pair_identity_verified.csv"

#: Which rows are displayed panels, and which are historical comparisons.
DISPLAYED_PANELS = {"class1_A", "class1_B", "class1_C", "class1_D", "class1_E",
                    "class1_F", "class2_A", "class2_B", "class2_C"}


def _s4b_frame() -> pd.DataFrame:
    """The recomputed superpositions, ordered as the approved table orders them."""
    from .structural import HISTORICAL_S4B_COLUMNS
    pub = _published("S4b")
    required = [c for c in pub.columns if c not in HISTORICAL_S4B_COLUMNS]
    d = _use_recomputed(TABLES / "structural_superpositions_recomputed.csv", required)
    if d is None:
        return None
    order = {n: i for i, n in enumerate(pub.panel_name)}
    return d.sort_values("panel_name", key=lambda s: s.map(order)).reset_index(drop=True)


def table_s4b() -> pd.DataFrame:
    """Per-pair superpositions under three named methods, plus domain RMSDs.

    **Current publication values only.** Fields describing the *submitted*
    version of the manuscript -- the allele its legend named, the RMSD it
    printed, and the deltas against it -- are not here; they are in
    ``TableS4b_provenance_submitted_version.csv``. Keeping them in this table
    invited exactly one misreading: that Figure 5F prints 3.4 A for HLA-C*07:01.
    It does not. The final figure prints 2.6 A for HLA-C*07:02, which is
    ``round(rmsd_cealign_CA, 1)`` of the recomputed measurement.

    ``record_level`` distinguishes the nine displayed panels from the one
    historical comparison row (``class1_F_alt``), which is a superposition
    against the allele the submitted legend named and is not a published panel.
    """
    from .structural import HISTORICAL_S4B_COLUMNS, IDENTITY_S4B_COLUMNS
    pub = _published("S4b")
    d = _s4b_frame()
    if d is None:
        # No PyMOL: the measurements come from the approved table, but the two
        # identity columns are taken from the verified alignment record rather
        # than from the approved table, which repeats class1_F's value on the
        # class1_F_alt row. Carrying the approved value here would reintroduce
        # the very error this release corrects.
        d = pub.copy()
        ver = VERIFIED_IDENTITY
        if ver.exists():
            v = pd.read_csv(ver).set_index("panel_name")
            for c in IDENTITY_S4B_COLUMNS:
                d[c] = d.panel_name.map(v[c])
            d["identity_method"] = ("carried from data/carried/"
                                    "tableS4b_pair_identity_verified.csv (MAFFT)")
    d = d.copy()
    d["record_level"] = d.panel_name.map(
        lambda n: "displayed panel" if n in DISPLAYED_PANELS
        else "historical comparison (not a displayed panel)")
    if "rmsd_displayed_in_final_figure" not in d.columns:
        d["rmsd_displayed_in_final_figure"] = d.rmsd_cealign_CA.round(1)
    d.loc[~d.panel_name.isin(DISPLAYED_PANELS), "rmsd_displayed_in_final_figure"] = pd.NA
    d["displayed_equals_round1_of_cealign"] = [
        (pd.notna(v) and v == round(float(c), 1)) if n in DISPLAYED_PANELS else pd.NA
        for n, v, c in zip(d.panel_name, d.rmsd_displayed_in_final_figure, d.rmsd_cealign_CA)]
    front = ["record_level", "final_figure", "panel", "panel_name", "mhc_class",
             "hla_allele_in_model", "sla_allele_in_model",
             "rmsd_displayed_in_final_figure", "displayed_equals_round1_of_cealign"]
    rest = [c for c in d.columns if c not in front and c not in HISTORICAL_S4B_COLUMNS]
    return d[front + rest]


def table_s4b_provenance() -> pd.DataFrame:
    """Historical record: what the SUBMITTED version of the manuscript carried.

    None of these values is a current publication value. They are retained
    because they document how the figures changed between submission and
    publication, and because Table S4b's ``delta_*`` columns are only
    interpretable beside them.
    """
    from .structural import HISTORICAL_S4B_COLUMNS
    pub = _published("S4b")
    d = _s4b_frame()
    src = pub if d is None else d.merge(
        pub[["panel_name", *[c for c in HISTORICAL_S4B_COLUMNS if c in pub.columns]]],
        on="panel_name", how="left", suffixes=("", "_pub"))
    out = pd.DataFrame({
        "panel_name": src.panel_name,
        "final_figure": src.final_figure,
        "panel": src.panel,
        "record_level": src.panel_name.map(
            lambda n: "displayed panel" if n in DISPLAYED_PANELS
            else "historical comparison (not a displayed panel)"),
        "hla_allele_in_model": src.hla_allele_in_model,
        "hla_allele_in_submitted_legend": src.get("hla_allele_in_legend"),
        "allele_match_vs_submitted_legend": src.get("allele_match"),
        "rmsd_printed_in_submitted_figure": src.get("rmsd_reported_in_manuscript"),
        "delta_cealign_minus_submitted": src.get("delta_cealign_minus_reported"),
        "delta_super_minus_submitted": src.get("delta_super_minus_reported"),
        "delta_align_minus_submitted": src.get("delta_align_minus_reported"),
    })
    out["note"] = ("Submitted-version record. The final manuscript and artwork print "
                   "round(rmsd_cealign_CA, 1); see Table S4b and "
                   "outputs/reports/structural_panel_audit.csv.")
    out.loc[out.panel_name == "class1_F_alt", "note"] = (
        "Historical comparison only. Superposition against HLA-C*07:01, the allele the "
        "SUBMITTED legend named. The final legend and artwork name HLA-C*07:02, which is "
        "the modelled allele in displayed panel class1_F. Not a published panel.")
    return out


def table_s5() -> pd.DataFrame:
    """Displayed structural positions with their annotations.

    A deterministic projection of ``structural_labels_v6.csv``: the displayed
    label set, with ``published_position`` presented as
    ``canonical_residue_number``. Labels use SLA/HLA residue order.
    """
    lab = pd.read_csv(SOURCE_DATA / "structural_labels_v6.csv")
    lab = lab.rename(columns={"published_position": "canonical_residue_number"})
    cols = list(_published("S5").columns)
    missing = [c for c in cols if c not in lab.columns]
    if missing:
        raise KeyError(f"Table S5 columns absent from structural_labels_v6.csv: {missing}")
    return lab[cols]


#: Pair counts the published coverage denominators use, per class.
S6_PAIR_DENOMINATOR = {"class1": 259, "class2": 151}


def table_s6() -> pd.DataFrame:
    """Experimental validation priorities. Mixed provenance, corrected coverage.

    The agreement half -- ``mean_agreement``, ``divergence``,
    ``n_pairs_covering`` and ``coverage`` -- is recomputed here from
    ``identity_eplet_agreement.parquet``. The eplet-name half comes from the
    licensed HLA Eplet Registry and is carried, as are the structural contact
    distances.

    **Coverage correction.** The approved table divides ``n_pairs_covering`` by
    the number of matrix rows present at that position rather than by the class
    pair count. Those coincide at every position covered by the full cohort, so
    526 of 527 rows are unaffected; class I position 312 appears in only the 22
    HLA-A pair rows and its coverage was 1/22 = 0.045455 instead of
    1/259 = 0.003861. The denominator is now fixed per class, so the
    transformation cannot reintroduce the error at another position.
    """
    pub = _published("S6")
    agr = pd.read_parquet(EPLET_AGREEMENT / "identity_eplet_agreement.parquet")
    frames = []
    for cls, denom in S6_PAIR_DENOMINATOR.items():
        a = agr[agr.analysis == cls]
        n_pairs = a.groupby(["sla", "hla"]).ngroups
        if n_pairs != denom:
            raise AssertionError(
                f"{cls}: the agreement matrix holds {n_pairs} pairs but the documented "
                f"coverage denominator is {denom}")
        g = a.groupby("position").agree_gap_skip
        frames.append(pd.DataFrame({
            "analysis": cls, "position": g.count().index,
            "mean_agreement_r": g.mean().values,
            "n_pairs_covering_r": g.count().values}))
    rec = pd.concat(frames, ignore_index=True)
    rec["divergence_r"] = 1 - rec.mean_agreement_r
    rec["coverage_r"] = rec.n_pairs_covering_r / rec.analysis.map(S6_PAIR_DENOMINATOR)

    out = pub.merge(rec, on=["analysis", "position"], how="left")
    missing = out.mean_agreement_r.isna()
    if missing.any():
        raise AssertionError(f"{int(missing.sum())} Table S6 positions absent from the "
                             "agreement matrix; cannot recompute them")
    for col in ("mean_agreement", "divergence", "n_pairs_covering", "coverage"):
        out[col] = out[f"{col}_r"]
    out["n_pairs_covering"] = out.n_pairs_covering.astype(int)
    return out[list(pub.columns)]


#: Historical S4b fields, renamed on the way into the provenance table.
_HIST_RENAME = {
    "hla_allele_in_legend": "hla_allele_in_submitted_legend",
    "allele_match": "allele_match_vs_submitted_legend",
    "rmsd_reported_in_manuscript": "rmsd_printed_in_submitted_figure",
    "delta_cealign_minus_reported": "delta_cealign_minus_submitted",
    "delta_super_minus_reported": "delta_super_minus_submitted",
    "delta_align_minus_reported": "delta_align_minus_submitted",
}

_HIST_DESC = {
    "hla_allele_in_submitted_legend":
        "HLA allele named in the SUBMITTED manuscript legend. HISTORICAL: the final "
        "legend names hla_allele_in_model.",
    "allele_match_vs_submitted_legend":
        "whether the model matches the SUBMITTED legend. HISTORICAL; not a statement "
        "about the final manuscript.",
    "rmsd_printed_in_submitted_figure":
        "RMSD printed in the SUBMITTED figure. HISTORICAL. The final figure prints "
        "round(rmsd_cealign_CA, 1); see Table S4b.",
    "delta_cealign_minus_submitted":
        "adopted CE value minus the value printed in the SUBMITTED figure. HISTORICAL.",
    "delta_super_minus_submitted":
        "adopted super value minus the value printed in the SUBMITTED figure. HISTORICAL.",
    "delta_align_minus_submitted":
        "adopted align value minus the value printed in the SUBMITTED figure. HISTORICAL.",
}

_S4B_NEW = [
    ("record_level", "object", '"displayed panel" (the nine published panels) or '
     '"historical comparison (not a displayed panel)" (class1_F_alt)'),
    ("rmsd_displayed_in_final_figure", "float64",
     "the value the FINAL panel prints: round(rmsd_cealign_CA, 1). Blank for the "
     "historical comparison row."),
    ("displayed_equals_round1_of_cealign", "bool",
     "check that the printed value equals the recomputed measurement rounded to one "
     "decimal; true for all nine displayed panels"),
    ("identity_method", "object",
     "how the two sequence-identity columns were produced in this run"),
]

_S4B_FIXED = {
    "percent_identity_over_shared_positions":
        "MAFFT identity between the two modelled sequences, over positions where both "
        "chains have a residue. RECOMPUTED here when MAFFT is present; otherwise carried "
        "from the approved table.",
    "n_positions_both_present":
        "alignment positions at which both modelled chains have a residue (MAFFT). "
        "RECOMPUTED here when MAFFT is present; otherwise carried.",
}


def data_dictionary() -> pd.DataFrame:
    """The approved dictionary, corrected so historical fields cannot be misread.

    Three changes, all confined to the structural tables:

    1. The six submitted-version fields move from the ``S4b`` block to a new
       ``S4b_provenance`` block and are renamed to say ``submitted``.
    2. Four new ``S4b`` entries describe what the final figure prints and how the
       identity columns were produced.
    3. The two identity descriptions name their method (MAFFT) and say whether
       this run recomputed them.

    The approved dictionary under ``publication_artwork/`` is not modified; the
    corrected one is written to ``outputs/``.
    """
    d = pd.read_csv(ARTWORK / "tables" / "PUBLICATION_TABLES_DATA_DICTIONARY.csv")
    hist = d[(d.table == "S4b") & (d.column.isin(_HIST_RENAME))].copy()
    d = d[~((d.table == "S4b") & (d.column.isin(_HIST_RENAME)))]
    hist["table"] = "S4b_provenance"
    hist["column"] = hist.column.map(_HIST_RENAME)
    hist["description"] = hist.column.map(_HIST_DESC)
    extra = pd.DataFrame([
        dict(table="S4b_provenance", column=c, dtype=t, n_non_null=10, example=e, description=desc)
        for c, t, e, desc in [
            ("panel_name", "object", "class1_F",
             "internal panel identifier; class1_F_alt is a historical comparison, not a panel"),
            ("final_figure", "object", "Figure 5", "final figure the panel belongs to"),
            ("panel", "object", "F", "panel letter"),
            ("record_level", "object", "displayed panel",
             '"displayed panel" or "historical comparison (not a displayed panel)"'),
            ("hla_allele_in_model", "object", "C*07:02",
             "HLA allele the model file actually contains"),
            ("note", "object", "Submitted-version record.",
             "why the row is retained and what it is not"),
        ]])
    new = pd.DataFrame([dict(table="S4b", column=c, dtype=t, n_non_null=10,
                             example="", description=desc) for c, t, desc in _S4B_NEW])
    d = d[~((d.table == "S4b") & (d.column.isin(set(new.column))))]
    for col, txt in _S4B_FIXED.items():
        d.loc[(d.table == "S4b") & (d.column == col), "description"] = txt
    out = pd.concat([d, new, hist, extra], ignore_index=True)
    order = {t: i for i, t in enumerate(
        ["S1", "S2", "S3", "S4a", "S4b", "S4b_provenance", "S5", "S6"])}
    return out.sort_values("table", key=lambda s: s.map(lambda t: order.get(t, 99)),
                           kind="stable").reset_index(drop=True)


def main(argv=None) -> int:
    ensure_outputs()
    out = TABLES / "publication"
    out.mkdir(parents=True, exist_ok=True)
    s1_published, s1_recomputed = table_s1()
    s1_published.to_csv(out / NAMES["S1"], index=False)
    s1_recomputed.to_csv(TABLES / "tableS1_haplotype_resolution_recomputed.csv", index=False)
    built = {"S2": table_s2(), "S3": table_s3(), "S4a": table_s4a(),
             "S4b": table_s4b(), "S4b_provenance": table_s4b_provenance(),
             "S5": table_s5(), "S6": table_s6()}
    for key, df in built.items():
        df.to_csv(out / NAMES[key], index=False)
    data_dictionary().to_csv(out / "PUBLICATION_TABLES_DATA_DICTIONARY.csv", index=False)
    eplets_recomputed = (TABLES / "eplet_statistics_v6.csv").exists()
    structural_recomputed = (TABLES / "structural_superpositions_recomputed.csv").exists() \
        and "rmsd_cealign_CA" in built["S4b"].columns \
        and (TABLES / "structural_model_confidence_recomputed.csv").exists()
    (TABLES / "table_provenance.json").write_text(json.dumps({
        "S1": "partial: haplotype columns regenerated from the curated workbooks; "
              "n_ciwd_hits/nearest_*/most_distant_* archived (need the upstream "
              "DIAMOND identity stage over IPD-IMGT/HLA and CIWD 3.0)",
        "S2": "recomputed",
        "S3": "recomputed" if eplets_recomputed else
              "carried from the archived statistics (the eplet position mask is not "
              "available locally; see slahla_pub.eplet_mask)",
        "S4a": "recomputed" if structural_recomputed else
               "archived (run `make structures` with PyMOL to recompute)",
        "S4b": "recomputed except percent_identity_over_shared_positions and "
               "n_positions_both_present, which are archived"
               if structural_recomputed else
               "archived (run `make structures` with PyMOL to recompute)",
        "S5": "recomputed", "S6": "archived (needs the licensed HLA Eplet Registry)",
    }, indent=2))
    print(f"tables: wrote {len(built) + 1} publication tables to {out}")
    print("  S1 partial, S2/S5 recomputed, S6 mixed (coverage corrected), S3 "
          + ("recomputed" if eplets_recomputed else "carried") + ", S4a/S4b "
          + ("recomputed from the deposited models"
             if structural_recomputed else "archived (no PyMOL; see `make structures`)"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
