#!/usr/bin/env python
"""Stage the archive deposit: build the tarballs, checksum them, nothing else.

This script does not upload. It writes `staging/` and a manifest; a human
uploads to the DOI archive and records the identifiers in PROVENANCE.md.

What goes in the deposit is what does not belong in a Git repository: the
recovered model checkpoints, and the vector artwork. Everything `make reproduce`
needs is already committed, so the deposit is for upstream reruns and for the
archival record the manuscript's availability statement promises.

    python scripts/stage_deposit.py --checkpoints /path/to/results/models

Archiving the Git repository does NOT archive the checkpoints: they are not in
it. Verify the deposited files list before citing a DOI.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import tarfile
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
STAGING = ROOT / "staging"
VERSION = "1.0.0"

#: The two recovered ORIGINAL seed-42 checkpoints. Seeds 43-46 are unavailable.
CHECKPOINTS = ["class1__armC_dedup_cluster", "class2__armC_dedup_cluster"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def add_tree(tar: tarfile.TarFile, src: Path, arc: str) -> int:
    n = 0
    for p in sorted(src.rglob("*")):
        if p.is_file():
            tar.add(p, arcname=f"{arc}/{p.relative_to(src)}")
            n += 1
    return n


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--checkpoints", type=Path,
                    help="directory holding the recovered checkpoint folders")
    ap.add_argument("--svg-source", type=Path,
                    help="directory of approved SVG artwork, if kept outside the repo")
    args = ap.parse_args()
    STAGING.mkdir(exist_ok=True)
    entries = []

    svg_src = args.svg_source or (ROOT / "outputs" / "figures")
    svgs = sorted(svg_src.glob("*.svg"))
    if svgs:
        p = STAGING / f"sla-hla-figures-svg-v{VERSION}.tar.gz"
        with tarfile.open(p, "w:gz") as t:
            for s in svgs:
                t.add(s, arcname=f"figures/{s.name}")
        entries.append({"file": p.name, "contents": f"{len(svgs)} vector figures",
                        "bytes": p.stat().st_size, "sha256": sha256(p)})

    if args.checkpoints:
        for name in CHECKPOINTS:
            d = args.checkpoints / name
            if not d.is_dir():
                print(f"  skipped (absent): {d}")
                continue
            p = STAGING / f"sla-hla-checkpoint-{name}-v{VERSION}.tar.gz"
            with tarfile.open(p, "w:gz") as t:
                n = add_tree(t, d, name)
            entries.append({
                "file": p.name,
                "contents": f"recovered ORIGINAL seed-42 checkpoint ({n} files)",
                "bytes": p.stat().st_size, "sha256": sha256(p),
                "provenance": "original training run, recovered; NOT a reconstruction"})
    else:
        print("  no --checkpoints given; checkpoint archives not staged")

    manifest = {
        "version": VERSION,
        "staged": date.today().isoformat(),
        "files": entries,
        "checkpoint_inventory": {
            "class1 seed 42": "ORIGINAL weights, recovered - included if staged",
            "class2 seed 42": "ORIGINAL weights, recovered - included if staged",
            "class1/class2 seeds 43-46": "ORIGINAL weights UNAVAILABLE; not reproducible "
                                         "by retraining (GPU kernels were not deterministic). "
                                         "The archived predictions from these runs ARE included "
                                         "in the repository and are what the published findings "
                                         "rest on.",
            "reconstructed checkpoints": "EXCLUDED. Later retraining produced different "
                                         "artefacts; they are not substitutes for the originals.",
        },
        "note": "Archiving the Git repository does not archive these files. Confirm the "
                "deposited file list before citing a DOI.",
    }
    (STAGING / "DEPOSIT_MANIFEST.json").write_text(json.dumps(manifest, indent=2))
    lines = [f"{e['sha256']}  {e['file']}" for e in entries]
    (STAGING / "SHA256SUMS.txt").write_text("\n".join(lines) + ("\n" if lines else ""))
    print(f"\nstaged {len(entries)} archive(s) in {STAGING}")
    for e in entries:
        print(f"  {e['file']}  {e['bytes'] / 1e6:.1f} MB")
    print(f"\nwrote {STAGING / 'DEPOSIT_MANIFEST.json'} and SHA256SUMS.txt")
    print("Nothing was uploaded.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
