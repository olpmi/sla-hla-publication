#!/usr/bin/env python
"""Package the candidate for independent review.

Contains the committed candidate (code, environments, packaged inputs,
documentation, approved reference artwork) plus the outputs of a full run, so a
reviewer can compare regenerated against approved without running anything first.

Excluded: the deposition staging directory and its checkpoint archives, Git
internals, caches, and anything from the private development repositories.

    python scripts/build_review_bundle.py
"""
from __future__ import annotations

import hashlib
import subprocess
import zipfile
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VERSION = "1.0.0"

#: ``data/local`` holds the user's own licensed registry derivative (the eplet
#: position mask). It is never bundled: see DATA_LICENSES.md.
EXCLUDE_DIRS = {".git", "staging", "__pycache__", ".pytest_cache", ".ruff_cache", "local"}
EXCLUDE_SUFFIX = {".pyc", ".pyo"}
EXCLUDE_NAMES = {".DS_Store"}

#: Outputs worth shipping. Everything else under outputs/ is intermediate.
OUTPUT_INCLUDE = (
    "outputs/figures",
    "outputs/tables/publication",
    "outputs/tables/SLA_Supplementary_Tables_corrected.xlsx",
    "outputs/reports",
)


def _git(*args: str) -> str:
    return subprocess.run(["git", "-C", str(ROOT), *args],
                          capture_output=True, text=True, check=True).stdout.strip()


def _keep(p: Path) -> bool:
    rel = p.relative_to(ROOT)
    if any(part in EXCLUDE_DIRS for part in rel.parts):
        return False
    if p.suffix in EXCLUDE_SUFFIX or p.name in EXCLUDE_NAMES:
        return False
    if rel.parts and rel.parts[0] == "outputs":
        return any(str(rel).startswith(k) for k in OUTPUT_INCLUDE)
    return True


def files() -> list[Path]:
    return sorted(p for p in ROOT.rglob("*") if p.is_file() and _keep(p))


def readme(commit: str, tag_commit: str, n_files: int) -> str:
    return f"""# Review bundle — SLA/HLA publication candidate {VERSION}

Built {date.today().isoformat()} from the local candidate. Nothing in this bundle
has been pushed, published or deposited.

    commit          {commit}
    tag v1.0.0  ->  {tag_commit}
    files           {n_files}

## What is where

| path | what it is |
|---|---|
| `README.md` | what the repository reproduces, how to install and run it |
| `REPRODUCIBILITY.md` | scientific starting points, upstream reruns, known limits |
| `VALIDATION_REPORT.md` | what was executed in both validation modes, what matched, what did not |
| `CORRECTION_LOG.md` | every value changed relative to the approved materials, with evidence |
| `PUBLICATION_OUTPUTS.csv` | every figure and table -> generator, inputs, verification |
| `MIRROR_INTO_SUBMISSION.md` | table and dictionary changes to apply to the delivered workbook |
| `DATA_LICENSES.md` | per-resource redistribution basis, including one that is unresolved |
| `src/`, `configs/`, `environment.yml`, `environments/` | the code and its environments |
| `data/` | packaged scientific inputs (read-only) |
| `tests/` | the retained test set |
| **`publication_artwork/`** | **the APPROVED reference figures and tables** — the comparison baseline, not regenerated |
| **`outputs/figures/`** | **the REGENERATED figures** from a full run of this candidate |
| **`outputs/tables/publication/`** | **the REGENERATED tables**, including the corrected data dictionary and `TableS4b_provenance_submitted_version.csv` |
| **`outputs/tables/SLA_Supplementary_Tables_corrected.xlsx`** | **the corrected workbook**; the delivered one under `publication_artwork/` is untouched |
| `outputs/reports/` | `structural_panel_audit.csv`, `tableS1_cell_differences.csv`, `tableS1_provenance.json`, `column_provenance.csv`, `missing_inputs.json`, `verification.json` |
| `data/carried/` | the explicit carried inputs: S4b alignment evidence, S1 haplotype selections |
| `SHA256SUMS.txt` | checksum of every file in this bundle |
| `COMMIT.txt` | the candidate commit and tag |

`publication_artwork/` and `outputs/` hold the same sixteen figures and the same
tables. The first is the approved artwork; the second is what this code produced.
They are kept apart deliberately so a reviewer can diff them.

## Reproducing from this bundle

```bash
micromamba env create -f environment.yml && micromamba activate slahla-pub
pip install -e .
make reproduce      # ~2-3 min, CPU only, no network
```

Optional, to recompute Tables S4a/S4b from the deposited models (needs PyMOL and
MAFFT):

```bash
micromamba env create -f environments/structures.yml
micromamba activate slahla-pub-structures
make structures && make reproduce
```

## The eplet position mask is not in this bundle

It is derived from the HLA Eplet Registry, whose terms prohibit redistributing
its tables and are silent on derivatives; no permission covering it has been
granted, so it is excluded rather than shipped on an unconfirmed reading. The
`outputs/figures/` in this bundle contain all sixteen figures because they were
produced with a locally supplied mask — running `make reproduce` from the bundle
alone produces thirteen, names the missing input and records it. See
`DATA_LICENSES.md` for the required registry release and input format, and
`VALIDATION_REPORT.md` for both validation modes.

## Not included

The deposition staging directory and its two ~1.5 GB checkpoint archives, Git
internals, caches, `data/local/` (any locally supplied licensed input), and every
private development and reviewer material. The checkpoint inventory and its
provenance are described in `PROVENANCE.md` and `REPRODUCIBILITY.md`.
"""


def main() -> int:
    commit = _git("rev-parse", "HEAD")
    tag_commit = _git("rev-list", "-n1", "v1.0.0")
    fs = files()
    stem = f"sla-hla-publication-v{VERSION}-review"
    out = ROOT.parent / f"{stem}.zip"

    lines, total = [], 0
    for p in fs:
        h = hashlib.sha256(p.read_bytes()).hexdigest()
        lines.append(f"{h}  {p.relative_to(ROOT)}")
        total += p.stat().st_size

    rd = readme(commit, tag_commit, len(fs) + 3)
    commit_txt = (f"candidate commit: {commit}\n"
                  f"tag v1.0.0 -> commit: {tag_commit}\n"
                  f"built: {date.today().isoformat()}\n"
                  f"files: {len(fs) + 3}\n"
                  f"pushed or published: no\n")

    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as z:
        for p in fs:
            z.write(p, f"{stem}/{p.relative_to(ROOT)}")
        z.writestr(f"{stem}/REVIEW_README.md", rd)
        z.writestr(f"{stem}/COMMIT.txt", commit_txt)
        z.writestr(f"{stem}/SHA256SUMS.txt", "\n".join(lines) + "\n")

    digest = hashlib.sha256(out.read_bytes()).hexdigest()
    print(f"{out}")
    print(f"  files      {len(fs) + 3}")
    print(f"  uncompressed {total / 1e6:.1f} MB")
    print(f"  zip size   {out.stat().st_size / 1e6:.1f} MB")
    print(f"  sha256     {digest}")
    print(f"  commit     {commit}")
    (ROOT.parent / f"{stem}.zip.sha256").write_text(f"{digest}  {out.name}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
