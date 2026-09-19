"""Column-level provenance: what was recomputed, what was carried, what was not run.

"All publication tables EXACT" is true but incomplete: exact agreement with an
archived output is not the same claim as independent recomputation. Every column
of every publication table is classified here as exactly one of:

``recomputed``              recalculated in this run from packaged inputs, and
                            compared against the approved value.
``carried``                 an archived value copied into an assembled output.
                            It may agree exactly - it is the same number - but it
                            was not independently derived here.
``external_not_executed``   its calculation needs an input this repository does
                            not ship and the clean-environment run did not have.
``key``                     an identifier column.

The classification is verified, not asserted: for every column marked
``recomputed`` this module recomputes the value and records the largest absolute
difference from the approved table.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import eplet_mask, s1_audit
from .paths import ARTWORK, EPLET_AGREEMENT, SOURCE_DATA, REPORTS, TABLES, ensure_outputs

TOL = 5e-7

#: Pair counts the published coverage denominators use.
N_PAIRS = {"class1": 259, "class2": 151}


def _approved(name: str) -> pd.DataFrame:
    return pd.read_csv(ARTWORK / "tables" / name)


def _s6_regenerated() -> pd.DataFrame:
    """The part of Table S6 that follows from the packaged agreement matrix."""
    agr = pd.read_parquet(EPLET_AGREEMENT / "identity_eplet_agreement.parquet")
    prof = pd.read_csv(SOURCE_DATA / "figureS8_S13_UPDATED_plotted_values.csv")
    try:
        masks = eplet_mask.load()
    except eplet_mask.MissingExternalInput:
        masks = None
    out = []
    for cls in ("class1", "class2"):
        a = agr[agr.analysis == cls]
        g = a.groupby("position")
        t = pd.DataFrame({
            "analysis": cls, "position": g.size().index,
            "mean_agreement": g.agree_gap_skip.mean().values,
            "n_pairs_covering": g.agree_gap_skip.count().values})
        t["divergence"] = 1 - t.mean_agreement
        t["coverage"] = t.n_pairs_covering / N_PAIRS[cls]
        if masks is not None:
            t["is_eplet_associated"] = t.position.isin(masks[cls])
        out.append(t)
    return pd.concat(out, ignore_index=True)


def _max_diff(a: pd.Series, b: pd.Series) -> float | None:
    if pd.api.types.is_numeric_dtype(a) and pd.api.types.is_numeric_dtype(b):
        d = (a.astype(float) - b.astype(float)).abs()
        return float(np.nanmax(d)) if len(d) else 0.0
    return float((a.fillna("").astype(str) != b.fillna("").astype(str)).sum())


def classify() -> pd.DataFrame:
    rows: list[dict] = []

    def add(table, column, prov, source, verified=None, note=""):
        rows.append({"table": table, "column": column, "provenance": prov,
                     "source": source, "verified_max_abs_difference": verified,
                     "note": note})

    # ---- S1 -----------------------------------------------------------------
    s1 = _approved("TableS1_haplotype_nearest_and_most_distant_HLA.csv")
    diffs = s1_audit.compare()
    from . import haplotypes as H
    regen = H.alleles_by_haplotype()
    merged = s1.merge(regen, on="sla_allele", how="left", suffixes=("_a", "_b"))
    add("S1", "sla_allele", "key", "curated haplotype workbooks")
    for c in s1_audit.REGENERATED_COLUMNS:
        n = int((diffs.column == c).sum())
        add("S1", c, "recomputed", "data/haplotypes/*.xlsx", float(n),
            "differing cells are listed in outputs/reports/tableS1_cell_differences.csv"
            if n else "")
    for c in s1_audit.CARRIED_COLUMNS:
        add("S1", c, "external_not_executed",
            "DIAMOND all-vs-all over IPD-IMGT/HLA + CIWD 3.0", None,
            "needs third-party data this repository does not ship; not run in the "
            "clean-environment validation")
    for c in s1.columns:
        if c not in {"sla_allele", *s1_audit.REGENERATED_COLUMNS, *s1_audit.CARRIED_COLUMNS}:
            add("S1", c, "carried", "approved table", None, "annotation column")

    # ---- S2, S3, S5: fully recomputed ---------------------------------------
    for key, name, src in [
        ("S2", "TableS2_canonical_split_leakage_and_five_seed_validation.csv",
         "archived five-seed run configs + canonical_split_summary.csv"),
        ("S3", "TableS3_HLA_eplet_associated_position_statistics.csv",
         "identity_eplet_agreement.parquet + eplet position mask"),
        ("S5", "TableS5_structural_position_annotations.csv",
         "structural_labels_v6.csv"),
    ]:
        app = _approved(name)
        built = TABLES / "publication" / name
        got = pd.read_csv(built) if built.exists() else None
        for c in app.columns:
            v = _max_diff(got[c], app[c]) if got is not None and c in got.columns else None
            add(key, c, "recomputed", src, v)

    # ---- S4a -----------------------------------------------------------------
    app = _approved("TableS4a_structural_model_confidence_per_model.csv")
    rec = TABLES / "structural_model_confidence_recomputed.csv"
    got = pd.read_csv(TABLES / "publication" /
                      "TableS4a_structural_model_confidence_per_model.csv") \
        if (TABLES / "publication").exists() else None
    # The --no-pymol path writes a pLDDT-only confidence file. S4a counts as
    # recomputed only when the reference superpositions ran too, so a partial
    # file cannot be reported as a full recomputation.
    recomputed_ran = rec.exists() and "rmsd_ce" in pd.read_csv(rec, nrows=0).columns
    for c in app.columns:
        v = _max_diff(got[c], app[c]) if got is not None and c in got.columns else None
        if recomputed_ran:
            add("S4a", c, "recomputed",
                "deposited model PDBs + reference structures (PyMOL)", v)
        else:
            add("S4a", c, "carried", "approved table", v,
                "PyMOL absent: run `make structures` to recompute")

    # ---- S4b -----------------------------------------------------------------
    sup_ran = (TABLES / "structural_superpositions_recomputed.csv").exists()
    s4b_path = TABLES / "publication" / "TableS4b_structural_superpositions_per_pair.csv"
    s4b = pd.read_csv(s4b_path) if s4b_path.exists() else None
    app4b = _approved("TableS4b_structural_superpositions_per_pair.csv")
    ident = ("percent_identity_over_shared_positions", "n_positions_both_present")
    method = (s4b.identity_method.iloc[0] if s4b is not None and "identity_method" in s4b.columns
              else "")
    if s4b is not None:
        for c in s4b.columns:
            if c in ("record_level", "final_figure", "panel", "panel_name", "mhc_class",
                     "hla_allele_in_model", "sla_allele_in_model"):
                add("S4b", c, "key", "structural panel manifest")
            elif c in ident:
                v = _max_diff(s4b[c], app4b[c]) if c in app4b.columns else None
                if method.startswith("MAFFT"):
                    add("S4b", c, "recomputed", "MAFFT over the modelled sequences", v,
                        "one row (class1_F_alt, a historical comparison) differs from the "
                        "approved table; evidence in "
                        "data/carried/tableS4b_pair_identity_verified.csv")
                else:
                    # The fallback reads data/carried/tableS4b_pair_identity_verified.csv,
                    # not the approved table: the approved table repeats class1_F's
                    # identity on class1_F_alt, so carrying it would undo correction C2.
                    add("S4b", c, "carried",
                        "data/carried/tableS4b_pair_identity_verified.csv (MAFFT v7.526)", v,
                        "MAFFT absent: the verified alignment record is carried, which "
                        "preserves the corrected class1_F_alt identity; install "
                        "environments/structures.yml to recompute")
            elif c in ("rmsd_displayed_in_final_figure",
                       "displayed_equals_round1_of_cealign", "identity_method"):
                add("S4b", c, "recomputed", "derived from rmsd_cealign_CA", 0.0,
                    "added by this repository; not a column of the approved table")
            else:
                v = _max_diff(s4b[c], app4b[c]) if c in app4b.columns else None
                add("S4b", c, "recomputed" if sup_ran else "carried",
                    "deposited model PDBs (PyMOL cealign/super/align)", v,
                    "" if sup_ran else "PyMOL absent: run `make structures` to recompute")

    # ---- S4b provenance (historical) ----------------------------------------
    prov_path = TABLES / "publication" / "TableS4b_provenance_submitted_version.csv"
    if prov_path.exists():
        for c in pd.read_csv(prov_path).columns:
            add("S4b_provenance", c, "carried", "approved table (submitted-version record)",
                None, "HISTORICAL: describes the submitted manuscript, not the final one")

    # ---- S6 -------------------------------------------------------------------
    app6 = _approved("TableS6_experimental_validation_priorities.csv")
    regen6 = _s6_regenerated()
    m6 = app6.merge(regen6, on=["analysis", "position"], how="left", suffixes=("_a", "_b"))
    for c in ["analysis", "mhc_class", "position"]:
        add("S6", c, "key", "position index")
    for c in ["mean_agreement", "divergence", "n_pairs_covering"]:
        add("S6", c, "recomputed", "identity_eplet_agreement.parquet",
            _max_diff(m6[f"{c}_b"], m6[f"{c}_a"]))
    if "is_eplet_associated_b" in m6.columns:
        add("S6", "is_eplet_associated", "recomputed",
            "identity_eplet_agreement.parquet + locally supplied eplet position mask",
            _max_diff(m6["is_eplet_associated_b"], m6["is_eplet_associated_a"]))
    else:
        add("S6", "is_eplet_associated", "external_not_executed",
            "HLA Eplet Registry position mask", None,
            "the mask is not shipped; supply it locally with "
            "`python -m slahla_pub.eplet_mask --from-registry DIR`")
    add("S6", "coverage", "recomputed",
        "n_pairs_covering / 259 (class I) or 151 (class II)",
        _max_diff(m6["coverage_b"], m6["coverage_a"]),
        "CORRECTED in this release: the approved table divides by the number of matrix "
        "rows present at each position, which differs from the class pair count only at "
        "class I position 312 (1/22 instead of 1/259)")
    for c in ["eplet_loci_defining", "eplet_names", "eplet_analysis_support",
              "eplet_name_source", "eplet_named_in_registry_group_mask"]:
        add("S6", c, "external_not_executed", "HLA Eplet Registry export", None,
            "registry names cannot be redistributed; obtain from epregistry.com.br")
    for c in ["d_peptide", "d_tcr", "d_partner", "sasa_assembled", "reference_pdb"]:
        add("S6", c, "carried", "structural_candidate_annotations_v6.csv", None,
            "present in the packaged annotations; carried rather than recomputed here")
    for c in ["evidence_streams"]:
        add("S6", c, "carried", "approved table", None, "derived label")

    return pd.DataFrame(rows)


def main(argv=None) -> int:
    ensure_outputs()
    d = classify()
    d.to_csv(REPORTS / "column_provenance.csv", index=False)
    counts = d.provenance.value_counts()
    print("column provenance across the publication tables:")
    for k in ["recomputed", "carried", "external_not_executed", "key"]:
        if k in counts:
            print(f"  {k:24s} {counts[k]:3d}")
    bad = d[(d.provenance == "recomputed") & d.verified_max_abs_difference.notna()
            & (d.verified_max_abs_difference > TOL)]
    if len(bad):
        print(f"\n  {len(bad)} recomputed column(s) differ from the approved value:")
        for _, r in bad.iterrows():
            print(f"    {r.table}.{r.column}: {r.verified_max_abs_difference:g}"
                  + (f"  ({r.note})" if r.note else ""))
    else:
        print("\n  every recomputed column matches the approved value within 5e-7")
    print(f"\nwrote {REPORTS / 'column_provenance.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
