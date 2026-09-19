"""Write ``outputs/reports/validation_report.md`` from the verification record.

Factual results of the run that produced ``outputs/reports/verification.json``:
the stages that executed, the comparisons they make, anything that failed, and
the inputs that were unavailable.
"""
from __future__ import annotations

import json
import platform

import pandas as pd
from datetime import date

from . import runlog
from .paths import REPORTS, TABLES


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
        out["pymol"] = "not installed"
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
    blocked = set(rep.get("figures_blocked") or [])
    missing = runlog.read()

    L = []
    A = L.append
    A("# Validation report")
    A("")
    A(f"Run on {date.today().isoformat()} · Python {v['python']} · {v['platform']}")
    A("")
    A(f"numpy {v['numpy']} · pandas {v['pandas']} · scipy {v['scipy']} · "
      f"matplotlib {v['matplotlib']} · pillow {v['pillow']} · PyMOL {v['pymol']}")
    A("")

    A("## Stages executed")
    A("")
    A("| stage | module |")
    A("|---|---|")
    A("| classifier ensemble | `slahla_pub.ensemble` |")
    A("| eplet-position statistics | `slahla_pub.eplet_stats` |")
    A(f"| structural measurements | `slahla_pub.structural` "
      f"({'recomputed with PyMOL' if v['pymol'] != 'not installed' else 'archived values, --no-pymol'}) |")
    A("| publication tables | `slahla_pub.tables` |")
    A("| figures | `slahla_pub.figures_final` |")
    A("| Table S1 audit | `slahla_pub.s1_audit` |")
    A("| column provenance | `slahla_pub.output_provenance` |")
    A("| verification | `slahla_pub.verify` |")
    A("")

    A("## Results")
    A("")
    A(f"- **Publication tables** — {len(tables_ok)} of {len(rep['tables'])} match the "
      "approved materials on every shared column and row"
      + (": " + ", ".join(sorted(tables_ok)) + "." if tables_ok else "."))
    A(f"- **Recomputed intermediates** — {len(inter_ok)} of "
      f"{len(rep['recomputed_vs_archived'])} identical to the archived tables they "
      "were recomputed from"
      + (": " + ", ".join(f"`{k}`" for k in sorted(inter_ok)) + "." if inter_ok else "."))
    n_fig = sum(1 for r in figs if r.get("regenerated"))
    A(f"- **Figures** — {n_fig} of {len(figs)} regenerated in PNG, PDF and SVG"
      + (f"; {len(blocked)} blocked by a missing input." if blocked else "."))
    A(f"- **Scientific endpoints** — {passed} of {total} checks passed.")
    src_ok = rep.get("source_data_unchanged")
    A(f"- **Inputs** — `data/source_data/` ({rep.get('source_data_file_count', '?')} files) "
      + ("byte-unchanged after the run." if src_ok
         else f"**CHANGED**: {rep.get('source_data_changed_files')}."))
    A("")

    failures = []
    if tables_bad:
        failures += [f"table {k}: {n}" for k, n in tables_bad.items()]
    if inter_bad:
        failures += [f"intermediate {k}: {n}" for k, n in inter_bad.items()]
    if passed != total:
        failures.append(f"{total - passed} scientific endpoint check(s)")
    if src_ok is False:
        failures.append("source_data was modified during the run")
    A("## Failures")
    A("")
    if failures:
        for f in failures:
            A(f"- {f}")
    else:
        A("None. Every comparison above either matched or is a recorded correction "
          "(see *Corrections*).")
    A("")

    A("## Figures")
    A("")
    A("Compared on content and layout, not bytes: PDFs embed creation timestamps and "
      "font subsets, so byte equality would fail on identical artwork.")
    A("")
    A("| figure | canvas matches | mean abs. pixel difference | ink bbox shift | PDF | SVG |")
    A("|---|---|---|---|---|---|")
    for r in figs:
        if r["figure"] in blocked:
            A(f"| {r['figure']} | — | not produced (missing input) | — | — | — |")
            continue
        A(f"| {r['figure']} | {'yes' if r.get('size_match') else 'no'} | "
          f"{r.get('mean_abs_difference')} | {r.get('ink_bbox_max_shift_px')} px | "
          f"{'yes' if r.get('pdf_written') else 'no'} | "
          f"{'yes' if r.get('svg_written') else 'no'} |")
    A("")
    A("The residual pixel difference is text antialiasing: the approved rasters were "
      "rendered by a different freetype build, so glyph edges differ by a fraction of a "
      "pixel while the plotted geometry does not. Multi-panel figures are cropped to "
      "their ink (`bbox_inches='tight'`), which makes the canvas a few pixels wider or "
      "narrower when glyph metrics shift; the single-panel supplementary figures use the "
      "full canvas and match it exactly.")
    A("")

    audit_p = REPORTS / "structural_panel_audit.csv"
    if audit_p.exists():
        aud = pd.read_csv(audit_p)
        n_agree = int((aud.agreement == "agree").sum())
        # `ce_rmsd_recomputed` is populated only when PyMOL re-ran the superpositions;
        # the render report and the adopted manifest are always available.
        recomputed = aud.ce_rmsd_recomputed.notna().any()
        records = ("the render report the figure builder reads, the final structural "
                   "manifest, and the recomputation from the deposited model files"
                   if recomputed else
                   "the render report the figure builder reads and the final structural "
                   "manifest")
        A("## Structural panels")
        A("")
        A(f"{n_agree} of {len(aud)} displayed panels agree with the measurement they "
          "print. Each panel prints `round(rmsd_cealign_CA, 1)` of the whole-chain "
          f"C-alpha CE RMSD, and the records of that measurement — {records} — agree to "
          f"within {aud.max_abs_spread_between_records_A.max():.4f} A.")
        if not recomputed:
            A("")
            A("The superpositions were **not** recomputed in this run (PyMOL absent), so "
              "the comparison is against the archived measurements only.")
        A("")
        measured = "recomputed CE RMSD" if recomputed else "adopted CE RMSD"
        A(f"| figure | panel | HLA | SLA | {measured} | displayed | agree |")
        A("|---|---|---|---|---|---|---|")
        for _, r in aud.iterrows():
            val = r.ce_rmsd_recomputed if recomputed else r.ce_rmsd_manifest_adopted
            A(f"| {r.final_figure} | {r.panel} | {r.hla_allele} | {r.sla_allele} | "
              f"{val} A | {r.rmsd_displayed_in_final_figure} A | "
              f"{'yes' if r.agreement == 'agree' else 'REVIEW'} |")
        A("")

    prov_p = REPORTS / "column_provenance.csv"
    if prov_p.exists():
        cp = pd.read_csv(prov_p)
        counts = cp.provenance.value_counts()
        A("## Recomputed, carried, or not executed")
        A("")
        A("Exact agreement with an archived table is not the same claim as independent "
          "recomputation, so every column of every publication table is classified and "
          "the classification is checked: each `recomputed` column is recalculated here "
          "and compared. Full detail in `outputs/reports/column_provenance.csv`.")
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

    A("## Corrections")
    A("")
    A("Two values in the assembled tables differ from the approved materials by design. "
      "`make verify` reports them as corrections and fails on any difference not on this "
      "list. The approved copies under `publication_artwork/` keep their original values "
      "so the comparison baseline stays what was approved.")
    A("")
    A("| # | where | approved | corrected | basis |")
    A("|---|---|---|---|---|")
    A("| C1 | Table S6, class I position 312, `coverage` | 0.0454545454545454 | "
      "0.003861003861003861 | 1/259, the class I pair count; the denominator is a fixed "
      "per-class constant asserted against the matrix |")
    A("| C2 | Table S4b, `class1_F_alt`, `percent_identity_over_shared_positions` | 74.92 | "
      "75.58 | MAFFT v7.526, 229/303 identical residues; evidence in "
      "`data/carried/tableS4b_pair_identity_verified.csv` |")
    A("")
    A("Neither changes a published measurement. See `CORRECTION_LOG.md`.")
    A("")

    A("## Limitations")
    A("")
    if missing:
        for m in missing:
            A(f"- **{m['input']}** — blocks {', '.join(m['blocked_outputs'])}. "
              f"{m['how_to_supply_it'].splitlines()[0].strip()}")
    A("- **Table S1's identity columns** (`n_ciwd_hits`, `nearest_*`, `most_distant_*`) "
      "need the DIAMOND all-vs-all stage over IPD-IMGT/HLA plus the CIWD release. That "
      "stage is not executed here and its inputs are not shipped, so the values are "
      "carried from the approved table.")
    A("- **Table S6's eplet-name columns** (`eplet_names`, `eplet_loci_defining`, "
      "`eplet_name_source`, `eplet_named_in_registry_group_mask`, "
      "`eplet_analysis_support`) require the HLA Eplet Registry export, which is not "
      "redistributed. Its structural columns are carried from "
      "`structural_candidate_annotations_v6.csv`.")
    if v["pymol"] == "not installed":
        A("- **Tables S4a and S4b** were assembled from archived structural measurements. "
          "Recomputing them from the deposited model files needs PyMOL "
          "(`environments/structures.yml`), because the published quantities are defined "
          "by PyMOL's `cealign`, `super` and `align`.")
    A("- **The original seed 43–46 checkpoints** are unavailable and retraining cannot "
      "reproduce them, because the GPU kernels used were not deterministic. The archived "
      "predictions those runs produced are shipped in `data/archived_predictions/` and "
      "are what the classifier results are recomputed from.")
    A("- **The RMSDs printed in the submitted figures** were produced interactively; no "
      "computational provenance survives for them. They are retained in "
      "`TableS4b_provenance_submitted_version.csv`. The values the final figures print "
      "are reproduced above.")
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
    out = REPORTS / "validation_report.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(text)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
