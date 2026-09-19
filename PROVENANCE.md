# Provenance

Where the code and the archived scientific inputs in this repository come from.

## Source

| | |
|---|---|
| Authoritative inputs | the author-approved manuscript package of 18 September 2026 |
| Package SHA-256 | `7b07cd3cd518d461083301bf6ff73f6630da87a6bd586b3ffb6b3d110978a66d` |
| Package integrity | all 174 files verified against the package's own `SHA256SUMS.txt` before export |

## How each module was derived

`slahla_pub` is a subset of the analysis code that produced the published
results, not a copy of it. It is imported under its own name so it cannot shadow,
or be shadowed by, the development package (`slahla`) on one `PYTHONPATH`.

Paths below that are not `slahla_pub` modules — `code/…`, `structural_render_inputs/…`,
`12_scripts_env/…` — are files of the approved manuscript package, not of this
repository.

- `figures_final.py` — from `code/build_figures_final.py` in the approved
  package, the authoritative source for the published artwork. Every
  measurement, colour, size, panel letter and layout constant is preserved; no
  plotted quantity changed.
- `figure_helpers.py` — `project` and `place_labels`, reproduced verbatim from
  `structural_render_inputs/figures_structural_final.py` in the approved package.
  That script reached them by `ast.parse` + `exec` because importing it pulls in
  PyMOL; they are ordinary functions here.
- `structural.py` — a faithful port of
  `12_scripts_env/recompute_structural_rmsd.py`, the authoritative generator for
  Table S4b. Only the paths were made configurable, and the two MAFFT identity
  columns were added from the method its data dictionary states.
- `ensemble.py` — the estimator definitions from `slahla.ensemble`, with
  canonical membership enforced rather than assumed.
- `eplet_stats.py` — the statistics from `slahla.eplet_identity`, extended to the
  v6 definitions (per-locus scopes, the `(extreme + 1)/(shifts + 1)` estimator).
- `haplotypes.py` — the workbook parser from `slahla.haplotypes`, unchanged.

## Scope of retained assets

Only what the retained analyses or their documented upstream reruns need:

- 10 archived prediction tables — the canonical five seeds × two classes.
- 24 model PDBs — exactly the set Tables S4a and S4b reference.
- 3 experimental reference structures — exactly those the published comparisons
  use (1AQD, 1HHK, 1S9V; RCSB PDB archive data are CC0).

`data/MANIFEST.json` records the source URL, SHA-256, byte count, retrieval
timestamp and release for the 54 third-party files the upstream workflow needs.

## Release

| | |
|---|---|
| Version | 1.0.0 |
| Repository | https://github.com/olpmi/sla-hla-publication |
| Release | https://github.com/olpmi/sla-hla-publication/releases/tag/v1.0.0 |
| Release tag | `v1.0.0` (annotated) |
| Commit | `git rev-parse v1.0.0` |
| Release assets | `SLA_Supplementary_Tables_corrected.xlsx` — the corrected supplementary workbook |
| Archive DOI | *pending — not yet issued* |
| Model checkpoints | *not deposited* — see below |
| Licence | MIT for code; see `DATA_LICENSES.md` for data |

**Archive DOI.** None has been issued. DOI deposition is a separate step from
this GitHub publication and has not been carried out. Nothing here cites a DOI.

**Model checkpoints.** The two recovered original seed-42 checkpoints (class I and
class II) are **not in this repository and have not been deposited**. Archiving
this Git repository does not archive them. Seeds 43–46 are unavailable and are not
reproducible by retraining, because the GPU kernels used were not deterministic;
checkpoints from later retraining are excluded as they are not substitutes for
the originals. The published findings rest on the archived predictions from all
five seeds, which **are** included here under `data/archived_predictions/`.
