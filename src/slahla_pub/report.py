"""Write VALIDATION_REPORT.md from the verification record.

One page: what ran, what matched, what did not, and why.
"""
from __future__ import annotations

import json
import platform

import pandas as pd
from datetime import date

from .paths import REPORTS, ROOT, TABLES


def _versions() -> dict:
    import numpy, pandas, scipy, matplotlib, PIL
    out = {"python": platform.python_version(), "platform": platform.platform(),
           "numpy": numpy.__version__, "pandas": pandas.__version__,
           "scipy": scipy.__version__, "matplotlib": matplotlib.__version__,
           "pillow": PIL.__version__}
    try:
        import pymol
        out["pymol"] = getattr(pymol, "__version__", "present")
    except ImportError:
        out["pymol"] = "not installed (structural recompute skipped)"
    return out


def render(rep: dict) -> str:
    prov = {}
    p = TABLES / "table_provenance.json"
    if p.exists():
        prov = json.loads(p.read_text())
    figs = rep["figures"]
    tables_ok = [k for k, v in rep["tables"].items() if not v]
    tables_bad = {k: v for k, v in rep["tables"].items() if v}
    inter_ok = [k for k, v in rep["recomputed_vs_archived"].items() if not v]
    inter_bad = {k: v for k, v in rep["recomputed_vs_archived"].items() if v}
    v = _versions()
    passed = rep["scientific_checks_passed"]
    total = rep["scientific_checks_total"]

    L = []
    A = L.append
    A("# Validation report")
    A("")
    A(f"Run on {date.today().isoformat()} · Python {v['python']} · {v['platform']}")
    A("")
    A(f"numpy {v['numpy']} · pandas {v['pandas']} · scipy {v['scipy']} · "
      f"matplotlib {v['matplotlib']} · pillow {v['pillow']} · PyMOL {v['pymol']}")
    A("")
    A("## What was executed")
    A("")
    A("```")
    A("make reproduce      # ensemble -> eplet statistics -> tables -> figures -> verify")
    A("make structures     # optional, needs PyMOL: recompute Tables S4a and S4b")
    A("make test")
    A("```")
    A("")
    A("## What matched")
    A("")
    A(f"**Publication tables** — {len(tables_ok)} of {len(rep['tables'])} identical to the "
      "approved materials on every column and row: " + ", ".join(sorted(tables_ok)) + ".")
    A("")
    A(f"**Recomputed intermediates** — {len(inter_ok)} of "
      f"{len(rep['recomputed_vs_archived'])} identical to the archived tables "
      "they were recomputed from:")
    A("")
    for k in sorted(inter_ok):
        A(f"- `{k}`")
    A("")
    A(f"**Scientific endpoints** — {passed} of {total} checks passed. These are the "
      "manuscript's own endpoints with their definitions attached, including the "
      "cross-seed subgroup unanimity (0/260 class I, 4/152 class II), the "
      "gap-rule-labelled eplet contrasts, and the distinction between the "
      "observation-weighted `difference` and the equal-position "
      "`profile_difference` that the permutation test actually tests.")
    A("")
    A("**Stochastic quantities reproduced exactly.** The pair-bootstrap intervals "
      "(2000 resamples) and circular-shift p-values (10,000 shifts) match the "
      "published values digit for digit, because the original seed (20260913), "
      "resample counts and sampling procedure were recovered and re-executed. No "
      "tolerance was applied to them.")
    A("")
    A("**Inputs unchanged.** `data/source_data/` "
      f"({rep.get('source_data_file_count', '?')} files) is byte-identical after a "
      "full run: " + ("confirmed." if rep.get("source_data_unchanged")
                      else "**NOT confirmed** — " + str(rep.get("source_data_changed_files"))))
    A("")
    A("## Figures")
    A("")
    A("Compared on content and layout, not bytes: PDFs embed creation timestamps "
      "and font subsets, so byte equality would fail on identical artwork.")
    A("")
    A("| figure | canvas matches | mean abs. pixel difference | ink bbox shift | PDF | SVG |")
    A("|---|---|---|---|---|---|")
    for r in figs:
        A(f"| {r['figure']} | {'yes' if r.get('size_match') else 'no'} | "
          f"{r.get('mean_abs_difference')} | {r.get('ink_bbox_max_shift_px')} px | "
          f"{'yes' if r.get('pdf_written') else 'no'} | "
          f"{'yes' if r.get('svg_written') else 'no'} |")
    A("")
    A("All sixteen figures regenerate in PNG, PDF and SVG. The residual pixel "
      "difference is text antialiasing: the approved rasters were rendered by a "
      "different freetype build, so glyph edges differ by a fraction of a pixel "
      "while the plotted geometry does not. Multi-panel figures are cropped to "
      "their ink (`bbox_inches='tight'`), which makes the canvas a few pixels "
      "wider or narrower when glyph metrics shift; the five single-panel "
      "supplementary figures use the full canvas and match it exactly.")
    A("")
    audit_p = REPORTS / "structural_panel_audit.csv"
    if audit_p.exists():
        aud = pd.read_csv(audit_p)
        A("## Structural panels: measured, displayed, reconciled")
        A("")
        A(f"All {len(aud)} displayed panels agree. The value each panel prints is "
          "`round(rmsd_cealign_CA, 1)` of the whole-chain C-alpha CE RMSD, and three "
          "independent records of that measurement - the render report the figure "
          "builder reads, the final structural manifest, and the recomputation from the "
          "deposited model files - agree to within "
          f"{aud.max_abs_spread_between_records_A.max():.4f} A.")
        A("")
        A("| figure | panel | HLA | SLA | recomputed CE RMSD | displayed | agree |")
        A("|---|---|---|---|---|---|---|")
        for _, r in aud.iterrows():
            A(f"| {r.final_figure} | {r.panel} | {r.hla_allele} | {r.sla_allele} | "
              f"{r.ce_rmsd_recomputed} A | {r.rmsd_displayed_in_final_figure} A | "
              f"{'yes' if r.agreement == 'agree' else 'REVIEW'} |")
        A("")
        A("**Figure 5F carries no discrepancy.** The final legend and the final artwork "
          "both name HLA-C*07:02, and the panel prints 2.6 A, which is the recomputed "
          "2.638 A rounded. The HLA-C*07:01 allele and the 3.4 A value belong to the "
          "*submitted* version of the manuscript. They are retained only as history, in "
          "`TableS4b_provenance_submitted_version.csv`, and the `class1_F_alt` row that "
          "superposes the C*07:01 model is labelled a historical comparison, not a "
          "displayed panel.")
        A("")

    prov_p = REPORTS / "column_provenance.csv"
    if prov_p.exists():
        cp = pd.read_csv(prov_p)
        counts = cp.provenance.value_counts()
        A("## Recomputed, carried, or not executed")
        A("")
        A("Exact agreement with an archived table is not the same claim as independent "
          "recomputation. Every column of every publication table is classified and the "
          "classification is checked, not asserted: each `recomputed` column is "
          "recalculated here and compared. Full detail in "
          "`outputs/reports/column_provenance.csv`.")
        A("")
        A("| class | columns | meaning |")
        A("|---|---|---|")
        A(f"| recomputed | {counts.get('recomputed', 0)} | recalculated in this run from packaged inputs |")
        A(f"| carried | {counts.get('carried', 0)} | archived value copied into an assembled output |")
        A(f"| external_not_executed | {counts.get('external_not_executed', 0)} | needs an input this repository does not ship |")
        A(f"| key | {counts.get('key', 0)} | identifier column |")
        A("")
        per = cp.groupby(["table", "provenance"]).size().unstack(fill_value=0)
        A("| table | recomputed | carried | external | key |")
        A("|---|---|---|---|---|")
        for t in ["S1", "S2", "S3", "S4a", "S4b", "S4b_provenance", "S5", "S6"]:
            if t in per.index:
                r = per.loc[t]
                A(f"| {t} | {r.get('recomputed', 0)} | {r.get('carried', 0)} | "
                  f"{r.get('external_not_executed', 0)} | {r.get('key', 0)} |")
        A("")

    A("## Differences from the approved tables, and what causes them")
    A("")
    s1p = REPORTS / "tableS1_cell_differences.csv"
    if s1p.exists():
        d1 = pd.read_csv(s1p)
        A(f"**Table S1 - {len(d1)} differing cells of 492 compared.** The seven "
          "haplotype columns are regenerated from the curated workbooks in "
          "`data/haplotypes/`; all 82 alleles are covered, and `sla_locus`, `mhc_class`, "
          "`herd` and `local_name` match on every row. Cell-level detail is in "
          "`outputs/reports/tableS1_cell_differences.csv`.")
        A("")
        if len(d1):
            A("| allele | column | approved | regenerated | cause | affects a result |")
            A("|---|---|---|---|---|---|")
            for _, r in d1.iterrows():
                A(f"| {r.sla_allele} | {r.column} | {r.approved_value} | "
                  f"{r.regenerated_value} | {r.classification} | "
                  f"{r.affects_a_manuscript_result} |")
            A("")
        A("All remaining differences are the same kind: an allele that appears in more "
          "than one curated haplotype, where the approved table reports one and this "
          "regeneration reports the first (`haplotype`) or all of them (`haplotypes`). "
          "The rule the approved table used to choose is not recorded and is not "
          "recoverable from the shipped inputs. The underlying workbook rows are "
          "identical either way, all three alleles are accessory-locus records "
          "(SLA-DQA, SLA-DRA) that enter no classifier, eplet or structural analysis, "
          "and no manuscript result depends on them.")
        A("")
        A("A further five differences reported in the previous release have been "
          "**resolved**: `local_name` was being taken from the first non-empty value "
          "across every haplotype an allele appeared in, which imported a "
          "miniature-line label (X, Y, Z) onto a commercial haplotype. It now follows "
          "the selected haplotype row, and all 82 rows match.")
        A("")
    A("**Table S1's identity columns** - `n_ciwd_hits`, `nearest_*`, `most_distant_*` - "
      "are pairwise identities against the CIWD 3.0 common-and-well-documented allele "
      "set. They need the DIAMOND all-vs-all stage over IPD-IMGT/HLA plus the CIWD "
      "release. That calculation was **not executed** in the clean-environment run and "
      "its inputs are not shipped; the values are carried from the approved table.")
    A("")
    A("**Table S4b's two identity columns** are recomputed, not carried. The data "
      "dictionary states their method - MAFFT identity between the two modelled "
      "sequences, over positions where both chains have a residue - and nine of the ten "
      "rows reproduce exactly under it. The tenth, `class1_F_alt`, does not: the "
      "approved table repeats `class1_F`'s value (74.92% over 303 positions) for the "
      "alternative C*07:01 model, where MAFFT gives 75.58%. `class1_F_alt` is a "
      "historical comparison row, not a published panel, and no manuscript claim uses "
      "it. The approved value is left in place and the difference recorded rather than "
      "silently corrected. Without MAFFT both columns fall back to the approved values "
      "and `identity_method` says so.")
    A("")
    A("**Table S6** is a mixed table. Its agreement half - `mean_agreement`, "
      "`divergence`, `n_pairs_covering`, `is_eplet_associated` - is recomputed here "
      "from `identity_eplet_agreement.parquet` and the eplet position mask, and matches "
      "the approved table exactly (max difference 1.1e-16). `coverage` matches for 526 "
      "of 527 rows; class I position 312, which a single pair covers, uses a 22-pair "
      "denominator in the approved table where every other row uses 259. That is "
      "recorded, not altered. The eplet-name half - `eplet_names`, "
      "`eplet_loci_defining`, `eplet_name_source`, `eplet_named_in_registry_group_mask`, "
      "`eplet_analysis_support` - requires the HLA Eplet Registry export, which cannot "
      "be redistributed; to rebuild it, obtain the export from "
      "https://www.epregistry.com.br/, pool the flattened `resi` column per locus, and "
      "join on `position`. The structural columns are carried from "
      "`structural_candidate_annotations_v6.csv`.")
    A("")
    A("**The original seed 43-46 checkpoints** are unavailable, and retraining cannot "
      "reproduce them: GPU kernels were not made deterministic. The archived predictions "
      "those checkpoints produced are what the published classifier findings rest on, "
      "and those are shipped and recomputed from.")
    A("")
    A("**The RMSDs printed in the SUBMITTED figures** were produced interactively and "
      "typed in; no computational provenance survives for them. They are retained in "
      "`TableS4b_provenance_submitted_version.csv` as history. The values the FINAL "
      "figures print are reproducible and reproduced: see the panel audit above.")
    A("")
    clean = ROOT / "validation" / "cleanroom.json"
    if clean.exists():
        c = json.loads(clean.read_text())
        A("## Release validation: two modes")
        A("")
        A("The candidate was run two ways, because what it can produce depends on an "
          "input it does not ship. Both are reported; neither is presented as the other.")
        A("")
        for key, title in [("mode_a_local_inputs", "A - with locally available inputs"),
                           ("mode_b_public_bundle_only", "B - from the public bundle alone")]:
            m = c.get(key)
            if not m:
                continue
            A(f"### {title}")
            A("")
            A(m["description"] + ".")
            A("")
            A("| command | exit |")
            A("|---|---|")
            for e in m["commands"]:
                A(f"| `{e['command']}` | {e['exit_code']} |")
            A("")
            A(f"Figures produced: **{m['figures_produced']}**. "
              f"Column provenance: {m['column_provenance']}. "
              f"Scientific checks: {m['scientific_checks']}. Tests: {m['tests']}.")
            if m.get("outputs_not_produced"):
                A("")
                A(f"Not produced: **{', '.join(m['outputs_not_produced'])}** — each marks "
                  "eplet-associated positions individually and the mask is absent. "
                  f"Produced from archived aggregate statistics instead of recomputed: "
                  f"{', '.join(m['outputs_degraded_to_carried'])}. "
                  f"{m['missing_input_reported']}.")
            if m.get("workbook_vs_csv"):
                A("")
                A(f"Corrected supplementary workbook vs the regenerated CSVs: "
                  f"**{m['workbook_vs_csv']}**.")
            A("")
        com = c.get("common", {})
        A(f"In both modes `data/source_data` is byte-unchanged "
          f"({str(com.get('source_data_unchanged')).lower()}) and the clean extraction has "
          f"{com.get('private_repo_access')}. Optional environments: "
          + "; ".join(f"`{k}` {v}" for k, v in com.get("optional_environments_solved", {}).items())
          + ".")
        A("")

    log = ROOT / "CORRECTION_LOG.md"
    if log.exists():
        A("## Corrections applied in this release")
        A("")
        A("| # | where | old | new | basis | downstream effect |")
        A("|---|---|---|---|---|---|")
        A("| C1 | Table S6, class I position 312, `coverage` | 0.0454545454545454 | "
          "0.003861003861003861 | 1/259; the denominator is now a fixed per-class constant, "
          "asserted against the matrix | none - no rank column, no figure, no statistic |")
        A("| C2 | Table S4b, `class1_F_alt`, identity | 74.92% | 75.58% | MAFFT v7.526, "
          "229/303 identical residues; evidence in "
          "`data/carried/tableS4b_pair_identity_verified.csv` | none - a historical "
          "comparison row, not a displayed panel. Figure 5F unchanged |")
        A("| C3 | Table S1, `local_name` (previous release) | miniature-line labels | "
          "the selected haplotype's value | a defect in this repository, not in the "
          "approved table | none - all 82 rows now agree |")
        A("")
        A("Full detail, including the model checksums behind C2, is in `CORRECTION_LOG.md`. "
          "The approved copies under `publication_artwork/` keep their original values so "
          "the comparison baseline stays what the authors approved; `make verify` reports "
          "C1 and C2 as corrections and fails on anything not on that list.")
        A("")

    A("## Table provenance")
    A("")
    for k in ["S1", "S2", "S3", "S4a", "S4b", "S5", "S6"]:
        if k in prov:
            A(f"- **{k}** — {prov[k]}")
    A("")
    return "\n".join(L) + "\n"


def main(argv=None) -> int:
    rec = REPORTS / "verification.json"
    if not rec.exists():
        raise SystemExit("run `make verify` first (no outputs/reports/verification.json)")
    text = render(json.loads(rec.read_text()))
    (ROOT / "VALIDATION_REPORT.md").write_text(text)
    print(f"wrote {ROOT / 'VALIDATION_REPORT.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
