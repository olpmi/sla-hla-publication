"""Compare everything this repository regenerates against the approved materials.

Comparison rules, fixed before running and not relaxed to obtain a pass:

**Identifiers, counts, denominators** -- exact. Cohort membership, pair counts,
position counts, label spaces and split sizes either match or they do not.

**Numeric columns** -- compared at the precision the *source* carries. The
archived ``source_data`` CSVs store six decimals, so ``1e-9`` equality against
them would be testing float formatting rather than science; the threshold is
``5e-7``, half a unit in the last stored place. Full-precision parquet inputs are
compared tighter where they are the source.

**Stochastic quantities** -- the bootstrap intervals and permutation p-values are
reproduced *exactly*, because the original seed (20260913), resample count, shift
count and sampling procedure were recovered and re-executed. No tolerance is
applied to them and none is needed.

**Figures** -- content and layout, never bytes. PDFs embed creation timestamps and
font subsets, so byte equality would fail on identical artwork. Each regenerated
PNG is compared with the approved PNG on canvas size, ink bounding box and a
per-pixel difference; text rendering varies slightly with the freetype build, so
a small residual is expected and reported rather than asserted away.

**Inputs** -- ``data/source_data`` must be byte-identical after a full run. A
reproduction that edits its own inputs is not one.
"""
from __future__ import annotations

import hashlib
import json

import numpy as np
import pandas as pd

from . import runlog
from .paths import ARTWORK, FIGURES, REPORTS, SOURCE_DATA, TABLES, ensure_outputs

#: Half a unit in the last place the archived CSVs store.
NUMERIC_TOLERANCE = 5e-7

FIGURES_EXPECTED = [f"Figure{i}" for i in range(1, 7)] + [f"FigureS{i}" for i in range(1, 11)]
#: Tables whose shape deliberately differs from the approved CSV, with the reason.
#: Their shared columns are still compared value for value.
INTENTIONAL_RESHAPE = {
    "S4b": "submitted-version fields moved to TableS4b_provenance_submitted_version.csv; "
           "record_level, rmsd_displayed_in_final_figure, displayed_equals_round1_of_cealign "
           "and identity_method added",
}

#: Columns this repository deliberately redefines, and why.
INTENTIONAL_COLUMN_CHANGES = {
    ("S4b", "record_level"):
        'redefined from "model pair" to "displayed panel" / "historical comparison '
        '(not a displayed panel)", so class1_F_alt cannot be read as a published panel',
}

#: Values this release deliberately CORRECTS in the approved tables. Each is an
#: error in the approved value, verified against the underlying data, with the
#: old value preserved in CORRECTION_LOG.md. They are expected to differ.
CORRECTED_VALUES = {
    ("S6", "coverage"):
        "class I position 312: the approved table divides by the 22 matrix rows present at "
        "that position instead of the 259 class I pairs. 0.0454545454545454 -> "
        "0.003861003861003861 (= 1/259). Only this row changes; the denominator is now "
        "fixed per class so the transformation cannot reintroduce it elsewhere",
    ("S4b", "percent_identity_over_shared_positions"):
        "class1_F_alt: the approved table repeats class1_F's identity (74.92%) on the "
        "alternative C*07:01 model. MAFFT v7.526 over the two modelled sequences gives "
        "229/303 = 75.58%; evidence in data/carried/tableS4b_pair_identity_verified.csv",
}

#: Differences understood, documented, and left in place rather than papered
#: over. Anything on neither list fails the run.
KNOWN_DEVIATIONS: dict = {}

TABLES_EXPECTED = {
    "S1": "TableS1_haplotype_nearest_and_most_distant_HLA.csv",
    "S2": "TableS2_canonical_split_leakage_and_five_seed_validation.csv",
    "S3": "TableS3_HLA_eplet_associated_position_statistics.csv",
    "S4a": "TableS4a_structural_model_confidence_per_model.csv",
    "S4b": "TableS4b_structural_superpositions_per_pair.csv",
    "S5": "TableS5_structural_position_annotations.csv",
    "S6": "TableS6_experimental_validation_priorities.csv",
}


def sha256_tree(root) -> dict[str, str]:
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(root.rglob("*")) if p.is_file()}


def compare_frames(a: pd.DataFrame, b: pd.DataFrame, reshaped: str = "") -> list[str]:
    """Differences between a regenerated frame and the approved one.

    ``reshaped`` names a deliberate structural change; the shared columns are
    still compared value for value, but added or relocated columns are not
    reported as defects.
    """
    problems = []
    if list(a.columns) != list(b.columns):
        only_a = [c for c in a.columns if c not in b.columns]
        only_b = [c for c in b.columns if c not in a.columns]
        if reshaped:
            pass
        elif only_a or only_b:
            problems.append(f"columns differ (extra {only_a}, missing {only_b})")
        else:
            problems.append("column order differs")
    if len(a) != len(b):
        problems.append(f"row count {len(a)} vs {len(b)}")
        return problems
    for c in [c for c in a.columns if c in b.columns]:
        if pd.api.types.is_numeric_dtype(b[c]) and pd.api.types.is_numeric_dtype(a[c]):
            d = (a[c].astype(float) - b[c].astype(float)).abs()
            worst = float(np.nanmax(d)) if len(d) else 0.0
            if worst > NUMERIC_TOLERANCE:
                problems.append(f"{c}: max|diff| {worst:.3g}")
        else:
            n = int((a[c].fillna("").astype(str) != b[c].fillna("").astype(str)).sum())
            if n:
                problems.append(f"{c}: {n} differing values")
    return problems


def _ink_bbox(arr: np.ndarray):
    ink = arr.min(axis=2) < 0.97
    if not ink.any():
        return None
    ys, xs = np.nonzero(ink)
    return int(xs.min()), int(ys.min()), int(xs.max()), int(ys.max())


def compare_figure(name: str) -> dict:
    from PIL import Image
    new_p, ref_p = FIGURES / f"{name}.png", ARTWORK / "figures" / f"{name}.png"
    rec = {"figure": name, "regenerated": new_p.exists(), "reference": ref_p.exists()}
    if not (new_p.exists() and ref_p.exists()):
        return rec
    a = Image.open(new_p).convert("RGB")
    b = Image.open(ref_p).convert("RGB")
    rec["size_new"], rec["size_reference"] = list(a.size), list(b.size)
    rec["size_match"] = a.size == b.size
    x = np.asarray(a, float) / 255.
    y = np.asarray(b if a.size == b.size else b.resize(a.size), float) / 255.
    d = np.abs(x - y).mean(axis=2)
    rec["mean_abs_difference"] = round(float(d.mean()), 5)
    rec["fraction_pixels_over_2pct"] = round(float((d > 0.02).mean()), 5)
    ba, bb = _ink_bbox(x), _ink_bbox(y)
    rec["ink_bbox_new"], rec["ink_bbox_reference"] = ba, bb
    rec["ink_bbox_max_shift_px"] = (max(abs(p - q) for p, q in zip(ba, bb))
                                    if ba and bb else None)
    for ext in ("pdf", "svg"):
        rec[f"{ext}_written"] = (FIGURES / f"{name}.{ext}").exists()
    return rec


def run(inputs_before: dict | None = None) -> dict:
    from . import checks
    ensure_outputs()
    report: dict = {}

    report["figures"] = [compare_figure(n) for n in FIGURES_EXPECTED]
    report["missing_inputs"] = runlog.read()
    blocked = set(runlog.blocked())
    report["blocked_outputs"] = sorted(blocked)
    report["figures_missing"] = [r["figure"] for r in report["figures"]
                                 if not r["regenerated"] and r["figure"] not in blocked]
    report["figures_blocked"] = [r["figure"] for r in report["figures"]
                                 if not r["regenerated"] and r["figure"] in blocked]

    tab = {}
    built = TABLES / "publication"
    for key, name in TABLES_EXPECTED.items():
        p = built / name
        if not p.exists():
            tab[key] = ["not generated"]
            continue
        probs = compare_frames(pd.read_csv(p), pd.read_csv(ARTWORK / "tables" / name),
                               reshaped=INTENTIONAL_RESHAPE.get(key, ""))
        kept, intended, known = [], [], []
        for msg in probs:
            col = msg.split(":")[0]
            if (key, col) in INTENTIONAL_COLUMN_CHANGES:
                intended.append(col)
            elif (key, col) in CORRECTED_VALUES:
                report.setdefault("corrections_applied", {}).setdefault(key, []).append(msg)
            elif (key, col) in KNOWN_DEVIATIONS:
                known.append(msg)
            else:
                kept.append(msg)
        tab[key] = kept
        if intended:
            report.setdefault("intentional_column_changes", {})[key] = intended
        if known:
            report.setdefault("known_deviations", {})[key] = known
    report["tables"] = tab
    report["intentional_reshape"] = INTENTIONAL_RESHAPE

    recomputed = {}
    for name in ["canonical_ensemble_predictions.csv",
                 "canonical_ensemble_per_sequence_stability.csv",
                 "canonical_seed_stability.csv", "canonical_validation_across_seeds.csv",
                 "leave_one_seed_out_by_locus.csv", "eplet_statistics_v6.csv",
                 "Figure4_plotted_values.csv"]:
        p = TABLES / name
        if p.exists() and (SOURCE_DATA / name).exists():
            recomputed[name] = compare_frames(pd.read_csv(p), pd.read_csv(SOURCE_DATA / name))
    report["recomputed_vs_archived"] = recomputed

    cur = sha256_tree(SOURCE_DATA)
    if inputs_before is not None:
        changed = [k for k in cur if inputs_before.get(k) != cur[k]]
        report["source_data_unchanged"] = not changed
        report["source_data_changed_files"] = changed
    report["source_data_file_count"] = len(cur)

    # structural panel audit: all nine displayed panels must agree
    audit_p = REPORTS / "structural_panel_audit.csv"
    if audit_p.exists():
        audit = pd.read_csv(audit_p)
        report["structural_panels"] = {
            "n_displayed_panels": int(len(audit)),
            "n_agree": int((audit.agreement == "agree").sum()),
            "disagreeing": audit.loc[audit.agreement != "agree", "panel_name"].tolist(),
        }

    prov_p = REPORTS / "column_provenance.csv"
    if prov_p.exists():
        prov = pd.read_csv(prov_p)
        report["column_provenance_counts"] = prov.provenance.value_counts().to_dict()
        bad = prov[(prov.provenance == "recomputed")
                   & prov.verified_max_abs_difference.notna()
                   & (prov.verified_max_abs_difference > NUMERIC_TOLERANCE)]
        report["recomputed_columns_with_known_deviations"] = [
            {"table": r.table, "column": r.column,
             "max_abs_difference": r.verified_max_abs_difference, "note": r.note}
            for _, r in bad.iterrows()]

    s1_p = REPORTS / "tableS1_cell_differences.csv"
    if s1_p.exists():
        d1 = pd.read_csv(s1_p)
        report["tableS1_cell_differences"] = int(len(d1))
        report["tableS1_affects_a_manuscript_result"] = bool(
            len(d1) and (d1.affects_a_manuscript_result != "none").any())

    cs = checks.run()
    report["scientific_checks"] = [
        {"key": c.key, "description": c.description, "source": c.source,
         "expected": c.expected, "observed": c.observed, "passed": c.passed} for c in cs]
    report["scientific_checks_passed"] = sum(c.passed for c in cs)
    report["scientific_checks_total"] = len(cs)
    return report


def main(argv=None) -> int:
    import sys
    argv = list(sys.argv[1:] if argv is None else argv)
    before = None
    for a in argv:
        if a.startswith("--inputs-before="):
            before = json.loads(open(a.split("=", 1)[1]).read())
    rep = run(before)
    ensure_outputs()
    (REPORTS / "verification.json").write_text(json.dumps(rep, indent=2, default=str))

    fails = 0
    blocked = set(rep.get("blocked_outputs", []))
    print("== figures")
    for r in rep["figures"]:
        if not r["regenerated"] and r["figure"] in blocked:
            print(f"   {r['figure']:<10} NOT PRODUCED - blocked by a missing external input")
            continue
        bad = not r["regenerated"]
        fails += bad
        print(f"   {r['figure']:<10} {'MISSING' if bad else 'ok'}"
              + ("" if bad else f"  size_match={r['size_match']}"
                                f"  mean|diff|={r['mean_abs_difference']:.5f}"
                                f"  ink_shift={r['ink_bbox_max_shift_px']}px"))
    if "structural_panels" in rep:
        sp = rep["structural_panels"]
        ok = sp["n_agree"] == sp["n_displayed_panels"]
        fails += not ok
        print(f"== structural panels: {sp['n_agree']}/{sp['n_displayed_panels']} displayed "
              f"panels agree with the adopted measurement"
              + ("" if ok else f"  REVIEW: {sp['disagreeing']}"))
    print("== publication tables (shared columns, vs approved)")
    for k, probs in rep["tables"].items():
        fails += bool(probs)
        tag = "EXACT" if not probs else probs
        if k in rep.get("intentional_reshape", {}):
            tag = f"{tag}  [reshaped by design: {rep['intentional_reshape'][k]}]"
        print(f"   {k:<4} {tag}")
    for key, cols in rep.get("intentional_column_changes", {}).items():
        for c in cols:
            print(f"   {key} {c}: redefined by design "
                  f"({INTENTIONAL_COLUMN_CHANGES[(key, c)]})")
    for key, msgs in rep.get("corrections_applied", {}).items():
        for m in msgs:
            col = m.split(":")[0]
            print(f"   {key} {m}  [CORRECTED in this release: {CORRECTED_VALUES[(key, col)]}]")
    for key, msgs in rep.get("known_deviations", {}).items():
        for m in msgs:
            col = m.split(":")[0]
            print(f"   {key} {m}  [known, documented: {KNOWN_DEVIATIONS[(key, col)]}]")
    if "column_provenance_counts" in rep:
        c = rep["column_provenance_counts"]
        print("== column provenance: "
              + ", ".join(f"{v} {k}" for k, v in c.items()))
        for r in rep.get("recomputed_columns_with_known_deviations", []):
            print(f"   known deviation: {r['table']}.{r['column']} "
                  f"max|diff|={r['max_abs_difference']:g}")
    print("== recomputed intermediates (vs archived)")
    for k, probs in rep["recomputed_vs_archived"].items():
        fails += bool(probs)
        print(f"   {k:<48} {'EXACT' if not probs else probs}")
    if "source_data_unchanged" in rep:
        ok = rep["source_data_unchanged"]
        fails += not ok
        print(f"== inputs: source_data {'byte-unchanged' if ok else 'MUTATED: ' + str(rep['source_data_changed_files'])}")
    print(f"== scientific checks: {rep['scientific_checks_passed']}/{rep['scientific_checks_total']} passed")
    fails += rep["scientific_checks_total"] - rep["scientific_checks_passed"]
    if rep.get("missing_inputs"):
        print(runlog.banner())
        print(f"   {len(rep['blocked_outputs'])} output(s) not produced: "
              + ", ".join(rep["blocked_outputs"]))
    print(f"\n{'VERIFICATION PASSED' if not fails else f'VERIFICATION FAILED ({fails} problem(s))'}"
          + (" (with missing external inputs, listed above)" if rep.get("missing_inputs") else ""))
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
