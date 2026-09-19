"""Write ``outputs/reports/validation_report.md`` from the verification record.

Factual results of the run that produced ``outputs/reports/verification.json``:
the comparisons it makes, anything that failed, and the inputs that were
unavailable.

The repository does not record which stages a given invocation executed, so this
module never claims one ran. Where a quantity could have come from either an
archived measurement or a recomputation, the provenance the pipeline itself
recorded decides which is reported -- not whether an optional dependency happens
to be importable.
"""
from __future__ import annotations

import json
import platform

import pandas as pd
from datetime import date

from . import runlog
from .paths import REPORTS, TABLES

#: Modules `make reproduce` wires together, in order. This is the declared
#: workflow, not evidence that any of it ran in a particular invocation.
WORKFLOW = [
    ("classifier ensemble", "slahla_pub.ensemble"),
    ("eplet-position statistics", "slahla_pub.eplet_stats"),
    ("structural measurements", "slahla_pub.structural"),
    ("publication tables", "slahla_pub.tables"),
    ("figures", "slahla_pub.figures_final"),
    ("Table S1 audit", "slahla_pub.s1_audit"),
    ("column provenance", "slahla_pub.output_provenance"),
    ("verification", "slahla_pub.verify"),
]


def _versions() -> dict:
    import numpy, pandas, scipy, matplotlib, PIL
    return {"python": platform.python_version(), "platform": platform.platform(),
            "numpy": numpy.__version__, "pandas": pandas.__version__,
            "scipy": scipy.__version__, "matplotlib": matplotlib.__version__,
            "pillow": PIL.__version__}


def failures(rep: dict) -> list[str]:
    """Everything ``verify.main`` counts as a problem, by the same criteria.

    Kept in step with ``verify.main``: a figure that was expected and not
    produced (blocked ones excluded there and here), a displayed structural
    panel that disagrees with its adopted measurement, a table or recomputed
    intermediate with differences, mutated inputs, and failed endpoint checks.
    """
    out: list[str] = []
    for name in rep.get("figures_missing") or []:
        out.append(f"figure {name}: expected and not produced")
    sp = rep.get("structural_panels")
    if sp and sp.get("n_agree") != sp.get("n_displayed_panels"):
        bad = ", ".join(sp.get("disagreeing") or []) or "unnamed"
        out.append(f"structural panels: {sp.get('n_agree')} of "
                   f"{sp.get('n_displayed_panels')} agree; disagreeing: {bad}")
    for k, probs in (rep.get("tables") or {}).items():
        if probs:
            out.append(f"table {k}: {probs}")
    for k, probs in (rep.get("recomputed_vs_archived") or {}).items():
        if probs:
            out.append(f"intermediate {k}: {probs}")
    if rep.get("source_data_unchanged") is False:
        out.append("inputs: data/source_data was modified during the run "
                   f"({rep.get('source_data_changed_files')})")
    short = rep.get("scientific_checks_total", 0) - rep.get("scientific_checks_passed", 0)
    if short:
        out.append(f"scientific endpoints: {short} check(s) failed")
    return out


def _structural_provenance(prov: dict, aud: pd.DataFrame | None) -> tuple[bool, bool]:
    """(measurements are archived, a recomputed column carries values).

    ``make reproduce`` invokes ``structural --no-pymol``, so the default
    workflow never recomputes the superpositions however PyMOL is installed.
    Values in ``ce_rmsd_recomputed`` therefore come from a separate
    ``make structures`` run, not from the invocation being reported.
    """
    archived = str(prov.get("S4a", "")).startswith("archived")
    has_recomputed = bool(
        aud is not None
        and "ce_rmsd_recomputed" in aud.columns
        and aud.ce_rmsd_recomputed.notna().any()
    )
    return archived, has_recomputed


def render(rep: dict) -> str:
    prov = {}
    p = TABLES / "table_provenance.json"
    if p.exists():
        prov = json.loads(p.read_text())
    figs = rep["figures"]
    tables_ok = [k for k, v in rep["tables"].items() if not v]
    inter_ok = [k for k, v in rep["recomputed_vs_archived"].items() if not v]
    v = _versions()
    passed = rep["scientific_checks_passed"]
    total = rep["scientific_checks_total"]
    blocked = set(rep.get("figures_blocked") or [])
    missing = runlog.read()
    fails = failures(rep)

    audit_p = REPORTS / "structural_panel_audit.csv"
    aud = pd.read_csv(audit_p) if audit_p.exists() else None
    s4_archived, has_recomputed = _structural_provenance(prov, aud)

    L = []
    A = L.append
    A("# Validation report")
    A("")
    A(f"Run on {date.today().isoformat()} · Python {v['python']} · {v['platform']}")
    A("")
    A(f"numpy {v['numpy']} · pandas {v['pandas']} · scipy {v['scipy']} · "
      f"matplotlib {v['matplotlib']} · pillow {v['pillow']}")
    A("")

    A("## Workflow components")
    A("")
    A("The modules `make reproduce` wires together. This is the declared workflow;"
      " the repository does not record which of it a given invocation ran.")
    A("")
    A("| component | module |")
    A("|---|---|")
    for label, mod in WORKFLOW:
        A(f"| {label} | `{mod}` |")
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
    fig_note = ""
    if blocked:
        fig_note += f"; {len(blocked)} blocked by a missing input"
    if rep.get("figures_missing"):
        fig_note += f"; {len(rep['figures_missing'])} missing"
    A(f"- **Figures** — {n_fig} of {len(figs)} regenerated in PNG, PDF and SVG"
      + fig_note + ".")
    if aud is not None:
        sp = rep.get("structural_panels") or {}
        A(f"- **Structural panels** — {sp.get('n_agree', '?')} of "
          f"{sp.get('n_displayed_panels', '?')} displayed panels agree with the "
          "measurement they print.")
    A(f"- **Scientific endpoints** — {passed} of {total} checks passed.")
    src_ok = rep.get("source_data_unchanged")
    A(f"- **Inputs** — `data/source_data/` ({rep.get('source_data_file_count', '?')} files) "
      + ("byte-unchanged after the run." if src_ok
         else f"**CHANGED**: {rep.get('source_data_changed_files')}."))
    A("")

    A("## Failures")
    A("")
    if fails:
        for f in fails:
            A(f"- {f}")
    else:
        A("None.")
    A("")

    A("## Figures")
    A("")
    A("Compared on content and layout, not bytes: PDFs embed creation timestamps and "
      "font subsets, so byte equality would fail on identical artwork.")
    A("")
    A("| figure | canvas matches | mean abs. pixel difference | ink bbox shift | PDF | SVG |")
    A("|---|---|---|---|---|---|")
    produced = []
    for r in figs:
        if r["figure"] in blocked:
            A(f"| {r['figure']} | — | not produced (blocked by a missing input) | — | — | — |")
            continue
        if not r.get("regenerated"):
            A(f"| {r['figure']} | — | **NOT PRODUCED** | — | — | — |")
            continue
        produced.append(r)
        A(f"| {r['figure']} | {'yes' if r.get('size_match') else 'no'} | "
          f"{r.get('mean_abs_difference')} | {r.get('ink_bbox_max_shift_px')} px | "
          f"{'yes' if r.get('pdf_written') else 'no'} | "
          f"{'yes' if r.get('svg_written') else 'no'} |")
    A("")
    if produced:
        diffs = [r.get("mean_abs_difference") or 0 for r in produced]
        shifts = [r.get("ink_bbox_max_shift_px") or 0 for r in produced]
        A(f"Largest mean absolute pixel difference {max(diffs):.5f}; largest ink "
          f"bounding-box shift {max(shifts)} px.")
        # Only offer the antialiasing account where the measurements support it:
        # geometry that has moved by more than a few pixels is not a glyph-edge effect.
        if max(diffs) <= 0.05 and max(shifts) <= 8:
            A("")
            A("Differences of this size, with the ink bounding box shifting by only a "
              "few pixels, are consistent with text antialiasing rather than a change "
              "in plotted geometry: the approved rasters were rendered by a different "
              "freetype build. Multi-panel figures are cropped to their ink "
              "(`bbox_inches='tight'`), which moves the canvas edge when glyph metrics "
              "shift.")
        A("")

    if aud is not None:
        sp = rep.get("structural_panels") or {}
        A("## Structural panels")
        A("")
        A(f"{sp.get('n_agree', '?')} of {sp.get('n_displayed_panels', len(aud))} "
          "displayed panels agree with the measurement they print. Each panel prints "
          "`round(rmsd_cealign_CA, 1)` of the whole-chain C-alpha CE RMSD, and the "
          "available records of that measurement agree to within "
          f"{aud.max_abs_spread_between_records_A.max():.4f} A.")
        A("")
        if s4_archived:
            A("`make reproduce` invokes `slahla_pub.structural --no-pymol`, so the "
              "superpositions are **not** recomputed by the default workflow. The "
              "measurements below are the archived adopted values."
              + ("  A `ce_rmsd_recomputed` column is present from a separate "
                 "`make structures` run and is shown alongside; it was not produced by "
                 "the invocation this report covers." if has_recomputed else ""))
        else:
            A("Tables S4a/S4b were recomputed from the deposited model files "
              "(`make structures`).")
        A("")
        cols = "| figure | panel | HLA | SLA | adopted CE RMSD |"
        sep = "|---|---|---|---|---|"
        if has_recomputed:
            cols += " recomputed CE RMSD |"
            sep += "---|"
        cols += " displayed | agree |"
        sep += "---|---|"
        A(cols)
        A(sep)
        for _, r in aud.iterrows():
            row = (f"| {r.final_figure} | {r.panel} | {r.hla_allele} | {r.sla_allele} | "
                   f"{r.ce_rmsd_manifest_adopted} A |")
            if has_recomputed:
                row += f" {r.ce_rmsd_recomputed} A |"
            row += (f" {r.rmsd_displayed_in_final_figure} A | "
                    f"{'yes' if r.agreement == 'agree' else 'REVIEW'} |")
            A(row)
        A("")

    prov_p = REPORTS / "column_provenance.csv"
    if prov_p.exists():
        cp = pd.read_csv(prov_p)
        counts = cp.provenance.value_counts()
        A("## Recomputed, carried, or not executed")
        A("")
        A("Exact agreement with an archived table is not the same claim as independent "
          "recomputation, so every column of every publication table is classified and "
          "the classification is checked: each `recomputed` column is recalculated and "
          "compared. Full detail in `outputs/reports/column_provenance.csv`.")
        A("")
        A("| class | columns | meaning |")
        A("|---|---|---|")
        A(f"| recomputed | {counts.get('recomputed', 0)} | recalculated from packaged inputs |")
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
    A("Neither correction changes a figure or manuscript conclusion. See "
      "`CORRECTION_LOG.md`.")
    A("")

    A("## Limitations")
    A("")
    if missing:
        for m in missing:
            A(f"- **{m['input']}** — blocks {', '.join(m['blocked_outputs'])}. "
              f"{m['how_to_supply_it'].splitlines()[0].strip()}")
    A("- **Table S1's identity columns** (`n_ciwd_hits`, `nearest_*`, `most_distant_*`) "
      "need the DIAMOND all-vs-all stage over IPD-IMGT/HLA plus the CIWD release. That "
      "stage is not part of `make reproduce` and its inputs are not shipped, so the "
      "values are carried from the approved table.")
    A("- **Table S6's eplet-name columns** (`eplet_names`, `eplet_loci_defining`, "
      "`eplet_name_source`, `eplet_named_in_registry_group_mask`, "
      "`eplet_analysis_support`) require the HLA Eplet Registry export, which is not "
      "redistributed. Its structural columns are carried from "
      "`structural_candidate_annotations_v6.csv`.")
    if s4_archived:
        A("- **Tables S4a and S4b** are assembled from archived structural measurements. "
          "Recomputing them from the deposited model files needs PyMOL "
          "(`environments/structures.yml`) and the separate `make structures` target, "
          "because the published quantities are defined by PyMOL's `cealign`, `super` "
          "and `align`.")
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
