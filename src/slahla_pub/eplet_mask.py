"""The HLA eplet-associated position mask: a user-supplied local input.

Why it is not in this repository
--------------------------------
The mask is derived from the HLA Eplet Registry export, whose terms are research
use only and prohibit redistribution of its tables. Whether a derived set of
integer residue positions falls under that prohibition is a question we put to
ourselves and could not answer: no permission has been granted, and no term
explicitly covers derivatives. An unresolved basis is not a basis, so the mask is
**excluded from the public candidate** rather than shipped on an interpretation
nobody has confirmed.

What it is
----------
61 mature-protein positions for class I and 52 for class II. Nothing else - no
eplet names, no allele assignments, no evidence categories.

Generating it locally
---------------------
With an appropriately obtained registry export::

    python -m slahla_pub.eplet_mask --from-registry /path/to/export

writes ``data/local/eplet_positions.csv`` (``analysis,position``), which the
workflow then picks up automatically. ``data/local/`` is git-ignored; nothing
written there is redistributed by this repository.

Required registry release and input format
------------------------------------------
Release: the HLA Eplet Registry export as accessed **2026-08-05** from
https://www.epregistry.com.br/ (an account is required). Later releases give a
different mask and therefore different statistics - the published numbers are
pinned to that access date.

**Supported format: the registry's per-allele export**, one CSV per database,
with these columns::

    Name,Evidence,Exposition,Allele,resi,resn
    1C,,High,C*01:02,[[1]],[['C']]
    2D,B,Intermediate,DQA1*01:01,[[2]],[['D']]

``Name`` is the *eplet* name (``62GE``, ``1C``) and carries no locus information;
``Allele`` is the allele the eplet occurs on and is where the locus comes from.
Confusing the two is the one mistake that silently yields an empty mask, so the
importer reads the locus from ``Allele`` only, and rejects a file that has no
``Allele`` column rather than guessing.

Database files may be named for their database (``ABC``, ``DRB``, ``DQ``,
``DP``) or anything else; the name is used only to report what was read, never to
assign a locus. Two files are enough for the published mask: the ABC export and
the DR/DQ/DP export.

The selection rule, unchanged from the published analysis
---------------------------------------------------------
1. Keep eplets with ``Evidence`` in {A1, A2, B} **or** ``Exposition == 'High'``
   -- the disjunction the original analysis applies, not a conjunction.
2. Keep rows whose ``Allele`` is at one of the class's loci: **A, B, C** for
   class I, and **DRB1, DRB3, DRB4, DRB5, DQB1** for class II. DP and the alpha
   chains are excluded: the SLA cohort has no comparator at those loci.
3. Flatten the nested ``resi`` column to integer positions and take the union
   per class.
4. Intersect with the positions the agreement matrix covers, so the mask cannot
   name a position the analysis never scores.

Without the mask
----------------
Three published figures cannot be produced at all - **Figures S8, S9 and S10** -
because each marks eplet-associated positions individually and there is nothing
to mark. Two more, **Figure 4 and Table S3**, still appear but from the archived
aggregate statistics rather than recomputed ones: those carry group means and
counts, not positions, so they ship. The workflow reports both cases explicitly
and records them; it does not skip anything silently. Everything else - Figures
1, 2, 3, 5, 6, S1-S7 and Tables S1, S2, S4a, S4b, S5, S6 - is unaffected.
"""
from __future__ import annotations

import re

import pandas as pd

from .paths import DATA, EPLET_AGREEMENT

MASK_PATH = DATA / "local" / "eplet_positions.csv"

#: The loci each class's comparison covers, as they appear in the export's
#: ``Allele`` column. Class II is the beta chains the cohort has comparators for;
#: DP and the alpha chains are deliberately absent.
CLASS_LOCI = {"class1": ("A", "B", "C"),
              "class2": ("DRB1", "DRB3", "DRB4", "DRB5", "DQB1")}

#: The published eplet filter: a disjunction, not a conjunction.
EVIDENCE_KEEP = frozenset({"A1", "A2", "B"})
EXPOSITION_KEEP = "High"

#: Columns the per-allele export must have.
REQUIRED_COLUMNS = ("Name", "Evidence", "Exposition", "Allele", "resi")

#: Published position counts, as a check on a locally generated mask.
EXPECTED_N = {"class1": 61, "class2": 52}


class MissingExternalInput(RuntimeError):
    """Raised when an output needs an input this repository does not ship."""

    def __init__(self, what: str, blocks: list[str], how: str):
        self.what, self.blocks, self.how = what, blocks, how
        super().__init__(f"{what}\n  blocks: {', '.join(blocks)}\n  {how}")


#: Outputs that cannot be produced at all without the mask: each marks
#: eplet-associated positions individually, so there is nothing to draw.
BLOCKED_OUTPUTS = ["FigureS8", "FigureS9", "FigureS10"]

#: Outputs that still appear, but from archived aggregate statistics rather than
#: recomputed ones. These carry no position information, so they ship; without
#: the mask they simply cannot be re-derived, and the provenance report says so.
DEGRADED_TO_CARRIED = ["Figure4", "TableS3"]


def available() -> bool:
    return MASK_PATH.exists()


def load() -> dict[str, set[int]]:
    """{analysis: {positions}}, or raise with instructions."""
    if not MASK_PATH.exists():
        raise MissingExternalInput(
            what=f"the HLA eplet position mask is not present ({MASK_PATH})",
            blocks=BLOCKED_OUTPUTS,
            how=("generate it locally from an appropriately obtained registry export:\n"
                 "    python -m slahla_pub.eplet_mask --from-registry /path/to/export\n"
                 "  see DATA_LICENSES.md and the module docstring for the required "
                 "release and input format"))
    d = pd.read_csv(MASK_PATH)
    return {c: set(int(x) for x in g.position) for c, g in d.groupby("analysis")}


def positions(which: str) -> set[int]:
    return load()[which]


def locus_of(allele: str) -> str:
    """``C*01:02`` -> ``C``. The locus lives in the allele, never in the eplet name."""
    return str(allele).split("*", 1)[0].strip()


def flatten_resi(value) -> list[int]:
    """The nested ``resi`` cell (``[[1]]``, ``[[62, 63]]``) as integer positions."""
    import ast
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    try:
        parsed = ast.literal_eval(str(value))
    except (ValueError, SyntaxError):
        return [int(m) for m in re.findall(r"\d+", str(value))]
    out: list[int] = []

    def walk(node):
        if isinstance(node, (list, tuple, set)):
            for item in node:
                walk(item)
        elif isinstance(node, bool):
            return
        elif isinstance(node, (int, float)):
            out.append(int(node))
        elif isinstance(node, str):
            out.extend(int(m) for m in re.findall(r"\d+", node))

    walk(parsed)
    return out


def read_export(export_dir) -> pd.DataFrame:
    """Every per-allele table under ``export_dir``, validated before use.

    Raises rather than returning something half-usable: a mask built from a
    misread export is worse than no mask, because it looks like a result.
    """
    from pathlib import Path
    export = Path(export_dir)
    if not export.is_dir():
        raise SystemExit(f"not a directory: {export}")
    candidates = [p for p in sorted(export.rglob("*")) if p.suffix.lower() == ".csv"]
    if not candidates:
        raise SystemExit(
            f"no .csv files under {export}. The supported input is the registry's "
            "per-allele export; see `python -m slahla_pub.eplet_mask --help` and "
            "DATA_LICENSES.md for the schema.")
    frames, rejected = [], []
    for p in candidates:
        try:
            d = pd.read_csv(p)
        except Exception as exc:
            rejected.append(f"{p.name}: unreadable ({exc.__class__.__name__})")
            continue
        missing = [c for c in REQUIRED_COLUMNS if c not in d.columns]
        if missing:
            rejected.append(f"{p.name}: missing column(s) {', '.join(missing)}")
            continue
        if d.empty:
            rejected.append(f"{p.name}: no rows")
            continue
        d = d.copy()
        d["source_file"] = p.name
        frames.append(d)
    if not frames:
        raise SystemExit(
            "no usable per-allele registry table found.\n  "
            + "\n  ".join(rejected)
            + f"\n\nExpected columns: {', '.join(REQUIRED_COLUMNS)}.\n"
              "`Name` is the eplet name and carries no locus; the locus is read from "
              "`Allele`. See DATA_LICENSES.md for the pinned release and format.")
    for note in rejected:
        print(f"  skipped {note}")
    return pd.concat(frames, ignore_index=True)


def from_registry(export_dir) -> pd.DataFrame:
    """Build the mask from a user-supplied registry export directory."""
    raw = read_export(export_dir)
    kept = raw[raw["Evidence"].isin(EVIDENCE_KEEP)
               | (raw["Exposition"].astype(str).str.strip() == EXPOSITION_KEEP)].copy()
    if kept.empty:
        raise SystemExit(
            "no eplet passed the published filter (Evidence in {A1, A2, B} or "
            "Exposition == 'High'). Check that the export carries its Evidence and "
            "Exposition columns; an export without them cannot reproduce the mask.")
    kept["locus"] = kept["Allele"].map(locus_of)

    agr = pd.read_parquet(EPLET_AGREEMENT / "identity_eplet_agreement.parquet")
    covered = {c: set(g.position) for c, g in agr.groupby("analysis")}

    rows, summary = [], []
    for cls, loci in CLASS_LOCI.items():
        sel = kept[kept["locus"].isin(loci)]
        if sel.empty:
            raise SystemExit(
                f"the export contains no allele at the {cls} loci "
                f"({', '.join(loci)}). Loci seen: "
                f"{', '.join(sorted(kept['locus'].unique()))}. The locus is read from the "
                "Allele column; a file whose Allele column holds eplet names instead "
                "cannot produce a mask.")
        pos: set[int] = set()
        for v in sel["resi"]:
            pos.update(flatten_resi(v))
        if not pos:
            raise SystemExit(
                f"{cls}: no residue positions parsed from the resi column. Expected "
                "nested lists such as [[62, 63]].")
        keep = sorted(pos & covered.get(cls, set()))
        rows += [{"analysis": cls, "position": x} for x in keep]
        summary.append((cls, len(sel), len(pos), len(keep)))

    for cls, n_rows, n_pos, n_keep in summary:
        note = ""
        if n_keep != EXPECTED_N[cls]:
            note = (f"  <-- expected {EXPECTED_N[cls]} for the pinned 2026-08-05 release; "
                    "a different release will differ")
        print(f"  {cls}: {n_rows} filtered eplet rows -> {n_pos} positions -> "
              f"{n_keep} covered{note}")
    return pd.DataFrame(rows)


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="Build the eplet position mask locally.")
    ap.add_argument("--from-registry", metavar="DIR",
                    help="directory of the registry's per-allele export CSVs "
                         "(columns: " + ", ".join(REQUIRED_COLUMNS) + ")")
    ap.add_argument("--check", action="store_true", help="report whether the mask is present")
    a = ap.parse_args(argv)
    if a.check or not a.from_registry:
        if available():
            m = load()
            print(f"eplet mask present: {MASK_PATH}")
            for c, s in m.items():
                print(f"  {c}: {len(s)} positions (expected {EXPECTED_N.get(c)})")
            return 0
        print(f"eplet mask ABSENT ({MASK_PATH})")
        print(f"  blocks: {', '.join(BLOCKED_OUTPUTS)}")
        print("  build it with: python -m slahla_pub.eplet_mask --from-registry DIR")
        return 1
    d = from_registry(a.from_registry)
    if d.empty:
        raise SystemExit("refusing to write an empty mask")
    MASK_PATH.parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(MASK_PATH, index=False)
    print(f"wrote {MASK_PATH}: {len(d)} positions")
    for c, g in d.groupby("analysis"):
        print(f"  {c}: {len(g)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
