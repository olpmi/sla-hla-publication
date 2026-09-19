"""Cell-level comparison of the regenerated half of Table S1 against the approved table.

Table S1 is a mixed table. Seven columns describe which SLA allele sits in which
haplotype, herd and local line: those are regenerated here from the curated
workbooks in ``data/haplotypes/``. Seven more are pairwise sequence identities
against the CIWD 3.0 allele set: those need the upstream DIAMOND stage over
IPD-IMGT/HLA and the CIWD release, are not executed here, and are carried from
the approved table.

This module compares only the regenerated half, cell by cell, and classifies
every difference.
"""
from __future__ import annotations

import pandas as pd

from . import haplotypes as H
import json

from .paths import ARTWORK, REPORTS, ensure_outputs

#: Columns this repository regenerates from the curated workbooks.
REGENERATED_COLUMNS = ["sla_locus", "mhc_class", "haplotype", "haplotypes", "herd", "local_name"]

#: Columns carried from the approved table; they need the upstream identity stage.
CARRIED_COLUMNS = ["n_ciwd_hits", "nearest_hla", "nearest_identity", "nearest_category",
                   "most_distant_hla", "most_distant_identity", "most_distant_category"]


def _cause(allele: str, column: str, approved: str, regenerated: str, rows: pd.DataFrame) -> tuple[str, str, str]:
    """(classification, cause, manuscript impact) for one differing cell."""
    n_hap = rows.haplotype.nunique()
    if column in ("haplotype", "haplotypes") and n_hap > 1:
        return (
            "selection rule",
            f"{allele} appears in {n_hap} haplotypes in the curated workbooks "
            f"({', '.join(sorted(rows.haplotype.unique()))}). The approved table reports one; "
            "this regeneration reports the first for `haplotype` and all of them for "
            "`haplotypes`. The rule the approved table used to pick one is not recorded and "
            "is not recoverable from the shipped inputs. The underlying workbook rows are "
            "identical either way - no allele, locus or herd assignment differs.",
            "none",
        )
    return ("unclassified", "not attributable to a known cause; inspect before release", "review")


def compare() -> pd.DataFrame:
    approved = pd.read_csv(ARTWORK / "tables" /
                           "TableS1_haplotype_nearest_and_most_distant_HLA.csv")
    regenerated = H.alleles_by_haplotype()
    parsed = H.load()
    merged = approved.merge(regenerated, on="sla_allele", how="left",
                            suffixes=("_approved", "_regenerated"))
    rows = []
    for col in REGENERATED_COLUMNS:
        a = merged[f"{col}_approved"].fillna("").astype(str).str.strip()
        b = merged[f"{col}_regenerated"].fillna("").astype(str).str.strip()
        for i in merged.index[a != b]:
            allele = merged.sla_allele[i]
            cls, cause, impact = _cause(allele, col, a[i], b[i], parsed[parsed.allele == allele])
            rows.append({
                "sla_allele": allele, "sla_locus": merged.sla_locus_approved[i],
                "column": col,
                "approved_value": merged[f"{col}_approved"][i],
                "regenerated_value": merged[f"{col}_regenerated"][i],
                "classification": cls, "cause": cause,
                "affects_a_manuscript_result": impact,
            })
    cols = ["sla_allele", "sla_locus", "column", "approved_value", "regenerated_value",
            "classification", "cause", "affects_a_manuscript_result"]
    return pd.DataFrame(rows, columns=cols)


def coverage() -> dict:
    approved = pd.read_csv(ARTWORK / "tables" /
                           "TableS1_haplotype_nearest_and_most_distant_HLA.csv")
    regenerated = H.alleles_by_haplotype()
    diffs = compare()
    n_cells = len(approved) * len(REGENERATED_COLUMNS)
    carried = H.carried_selection()
    return {
        "approved_rows": int(len(approved)),
        "alleles_covered_by_the_workbooks": int(approved.sla_allele.isin(regenerated.sla_allele).sum()),
        "columns_verified_against_the_approved_table": REGENERATED_COLUMNS,
        "columns_requiring_an_external_input": CARRIED_COLUMNS,
        "external_input_required": "DIAMOND all-vs-all over IPD-IMGT/HLA and CIWD 3.0",
        "cells_compared": int(n_cells),
        "cells_differing": int(len(diffs)),
        "alleles_with_a_carried_haplotype_selection": sorted(carried.sla_allele.tolist()),
        "n_values_carried": int(2 * len(carried)),
        "any_difference_affects_a_manuscript_result":
            bool(len(diffs) and (diffs.affects_a_manuscript_result != "none").any()),
    }


def main(argv=None) -> int:
    ensure_outputs()
    d = compare()
    d.to_csv(REPORTS / "tableS1_cell_differences.csv", index=False)
    c = coverage()
    (REPORTS / "tableS1_provenance.json").write_text(json.dumps(c, indent=2))
    print(f"Table S1: {c['alleles_covered_by_the_workbooks']}/{c['approved_rows']} alleles "
          f"covered by the curated workbooks")
    print("  columns verified against the approved table: "
          + ", ".join(c["columns_verified_against_the_approved_table"]))
    print("  columns requiring an external input (NOT verified here): "
          + ", ".join(c["columns_requiring_an_external_input"]))
    print(f"  cells compared: {c['cells_compared']}, differing: {c['cells_differing']}")
    if c["alleles_with_a_carried_haplotype_selection"]:
        print(f"  haplotype selections carried from the approved table "
              f"({c['n_values_carried']} values over "
              f"{len(c['alleles_with_a_carried_haplotype_selection'])} alleles): "
              + ", ".join(c["alleles_with_a_carried_haplotype_selection"]))
    if len(d):
        for _, r in d.iterrows():
            print(f"  {r.sla_allele:18s} {r.column:11s} approved={r.approved_value!r} "
                  f"regenerated={r.regenerated_value!r} [{r.classification}]")
    print(f"  affects a manuscript result: "
          f"{c['any_difference_affects_a_manuscript_result']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
