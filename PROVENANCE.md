# Provenance

This repository is a minimal publication subset, exported from a private
development record that remains intact and unmodified.

## Source

| | |
|---|---|
| Authoritative inputs | `SLA_HLA_Manuscript.zip` → `SLA_Final_Manuscript_Package.zip`, the author-approved manuscript delivery of 18 September 2026 |
| Package SHA-256 | `7b07cd3cd518d461083301bf6ff73f6630da87a6bd586b3ffb6b3d110978a66d` |
| Package integrity | all 174 files verified against the package's own `SHA256SUMS.txt` before export |
| Development repository (revision) | `sla-hla-manuscript` at commit `56a154f6966f9b24027c4423da4a0370e5e56024` |
| Development repository (original) | `sla-hla-comparison` at commit `14253e02bf1e5b51b6ecb6edc0aaeb8886a01843` |
| Reviewed publication candidate | private working repository `sla-hla-publication` at commit `69233aa20275e99da19f82f43f1a997b37905f82` |

Both development repositories are private and were not modified, re-pointed or
re-published by this export. This repository has fresh Git history; no
development `.git` directory was copied.

## The public history

This public repository was exported from the reviewed candidate commit
`69233aa20275e99da19f82f43f1a997b37905f82` and opens with a **single initial
commit**. The candidate's own iterative correction history is deliberately not
published: earlier commits in it carry inputs that were subsequently removed from
the public candidate on the licensing grounds set out in `DATA_LICENSES.md`, and
publishing that history would redistribute them.

Every file here is byte-identical to that reviewed commit, with the exception of
this section and the repository URL and citation metadata added to `README.md`,
`CITATION.cff` and the release table below. No code, data, figure or result
differs. The reviewed candidate remains intact in its private working repository.

## Relationship to the development package

The development package is imported as `slahla`. This one is `slahla_pub`, so the
two cannot shadow each other on one `PYTHONPATH` and no result can silently come
from the wrong tree.

`slahla_pub` is a subset, not a copy. What was carried over:

- `figures_final.py` — from `code/build_figures_final.py` **in the delivery
  package** (not a file of this repository), which is the authoritative source
  for the published artwork. Five packaging repairs are documented in its module
  docstring; no plotted quantity changed.
- `figure_helpers.py` — `project` and `place_labels`, reproduced verbatim from
  `structural_render_inputs/figures_structural_final.py` **in the delivery
  package** (that module is not part of this repository). The delivery script
  reached them by `ast.parse` + `exec` because importing it pulls in PyMOL; they
  are ordinary functions here.
- `structural.py` — the superposition code is a faithful port of
  `12_scripts_env/recompute_structural_rmsd.py` **in the manuscript support
  package** (not a file of this repository), which is the authoritative generator
  for Table S4b. Only the paths were made configurable, and the two MAFFT
  identity columns were added from the method its data dictionary states.
- `ensemble.py` — the estimator definitions from `slahla.ensemble`, with
  canonical membership enforced rather than assumed.
- `eplet_stats.py` — the statistics from `slahla.eplet_identity`, extended to the
  v6 definitions (per-locus scopes, the `(extreme + 1)/(shifts + 1)` estimator).
- `haplotypes.py` — the workbook parser from `slahla.haplotypes`, unchanged.

## What was deliberately left behind

The integrated-gradients attribution pipeline and the residue-replacement
(perturbation) pipeline, with their figures, tables and tests: both analyses are
withdrawn from the publication. The Word tracked-changes machinery and the
manuscript-editing infrastructure. Handoff and packaging builders, transport
archives and review packages. Reviewer correspondence and earlier manuscript
drafts. Agent prompts, session transcripts, planning documents and debugging
reports. Development caches and generated PyMOL scripts carrying absolute
personal paths.

A test asserts the withdrawn modules are not importable, so the workflow cannot
come to depend on one by accident.

## Scope of retained assets

Only what the retained analyses or their documented upstream reruns need:

- 10 archived prediction tables (the canonical five seeds × two classes), not the
  16 present in the development tree.
- 24 model PDBs, not the 121 in the development tree — exactly the set Tables S4a
  and S4b reference.
- 3 experimental reference structures, not 11 — exactly those the published
  comparisons use.

## Release

| | |
|---|---|
| Version | 1.0.0 |
| Release notes | see `git tag -n99 v1.0.0` and VALIDATION_REPORT.md |
| Repository URL | https://github.com/olpmi/sla-hla-publication |
| Release | https://github.com/olpmi/sla-hla-publication/releases/tag/v1.0.0 |
| Release tag | `v1.0.0` (annotated, at the initial public commit) |
| Commit | `git rev-parse v1.0.0` |
| Release assets | `SLA_Supplementary_Tables_corrected.xlsx` — the corrected submission workbook |
| Archive DOI | *pending — not yet issued* |
| Model checkpoints | *not deposited* — see below |
| Licence | MIT for code; see `DATA_LICENSES.md` for data |

**Archive DOI.** None has been issued. DOI deposition is a separate step from
this GitHub publication and has not been carried out. Nothing in this repository
cites a DOI, and none has been invented.

**Model checkpoints.** The two recovered ORIGINAL seed-42 checkpoints (class I
and class II) are **not in this repository and have not been uploaded anywhere**.
Archiving this Git repository does not archive them. Seeds 43–46 are unavailable
and are not reproducible by retraining, because the GPU kernels used were not
deterministic; reconstructed checkpoints are excluded as they are not substitutes
for the originals. The published findings rest on the archived predictions from
all five seeds, which **are** included here under `data/archived_predictions/`.

The manuscript's Data and Code Availability statement can now carry the
repository URL, release tag and commit identifier from this table. The archive
DOI remains a placeholder until deposition is complete.
