# Comparative analysis of human HLA and swine SLA MHC sequences — reproduction package

Code, archived scientific inputs and documentation to regenerate every figure,
table and reported statistic of *"Comparative Bioinformatic Analysis of Human and
Swine MHC Sequences: Implications for Xenotransplantation"*.

**Repository:** https://github.com/olpmi/sla-hla-publication · **release:** `v1.0.0`

## What this reproduces

| | |
|---|---|
| **Main figures** | 1–6 |
| **Supplementary figures** | S1–S10 |
| **Supplementary tables** | S1, S2, S3, S4a, S4b, S5, S6, plus an S4b provenance sheet |
| **Statistics** | five-seed classifier summaries, eplet-position contrasts with their bootstrap intervals and permutation tests, structural superpositions and model confidence |

`PUBLICATION_OUTPUTS.csv` gives every output its generator, inputs, upstream
dependency and verification method. `outputs/reports/column_provenance.csv`
classifies every *column* of every table as independently recomputed, carried
from the approved materials, or dependent on an external input that was not
executed here — because exact agreement with an archived table is not the same
claim as independent recomputation.

## Installation

```bash
micromamba env create -f environment.yml     # or: conda env create -f environment.yml
micromamba activate slahla-pub
pip install -e .
```

Python ≥ 3.11. The core environment is CPU-only: numpy, pandas, scipy,
matplotlib, pillow, pyarrow, openpyxl, pyyaml.

## Required external input: the HLA eplet position mask

Supply this **before** reproducing, or three figures cannot be drawn.

The mask is derived from the HLA Eplet Registry, whose terms prohibit
redistributing its tables and say nothing about derivatives; no permission
covering this one has been granted, so it is not shipped. Obtain the registry's
per-allele CSV export from https://www.epregistry.com.br/ (account required,
release as accessed **2026-08-05**), then:

```bash
python -m slahla_pub.eplet_mask --from-registry /path/to/export
make mask        # confirm it was picked up
```

The export carries `Name,Evidence,Exposition,Allele,resi,resn`; the locus is read
from `Allele`, never from the eplet name. Against the pinned release the importer
reproduces the published mask exactly — the same 61 class I and 52 class II
positions, by membership. Anything written to `data/local/` stays there and is
not redistributed.

Without the mask the bundle produces **13 of the 16 figures**: Figures S8, S9 and
S10 cannot be drawn at all, and Figure 4 and Table S3 come from archived
aggregate statistics rather than recomputed ones. `make reproduce` prints a
`MISSING EXTERNAL INPUTS` block naming what is absent and what it blocks, and
records it in `outputs/reports/missing_inputs.json`.

Rerunning the optional upstream stages additionally needs IPD-IMGT/HLA and
IPD-MHC sequences and CIWD 3.0.0, which are likewise not redistributed.
`data/MANIFEST.json` pins the source URL, SHA-256, byte count, retrieval date and
release for every such file; `REPRODUCIBILITY.md` gives the acquisition
procedure. See `DATA_LICENSES.md` before redistributing anything in `data/`.

## Reproduce the published results

```bash
make reproduce
```

No GPU, no network, no retraining, no structure prediction. It recomputes the
classifier ensemble from the archived five-seed predictions, recomputes the eplet
statistics from the archived per-pair agreement matrix, assembles the publication
tables, regenerates the figures, then verifies everything against the approved
materials.

**Runtime** roughly 2–3 minutes on a laptop; ~1 GB RAM. The eplet permutation
test (10,000 circular shifts × 14 scope/reading combinations) is the slow part at
about 30 seconds.

**Outputs** land under `outputs/`:

```
outputs/figures/            Figure1-6, FigureS1-S10 in PNG, PDF and SVG
outputs/tables/publication/ TableS1, S2, S3, S4a, S4b, S5, S6
outputs/tables/             recomputed intermediates (ensemble, eplet statistics)
outputs/tables/SLA_Supplementary_Tables_corrected.xlsx   the corrected workbook
outputs/audit/              panel-letter validation, linkage orders, edge audits
outputs/reports/            verification.json, column_provenance.csv,
                            validation_report.md, missing_inputs.json
```

`data/` is read-only input and `make verify` asserts it is byte-unchanged after a
run. `make clean` removes `outputs/` and nothing else.

### The two supplementary workbooks

| file | what it is |
|---|---|
| `outputs/tables/SLA_Supplementary_Tables_corrected.xlsx` | **The corrected supplementary workbook — the one to use.** Built by `make workbook` (also run by `make reproduce`), it carries the two value corrections and the new `S4b_provenance` sheet. It is attached to the `v1.0.0` release so it can be had without running the pipeline. |
| `publication_artwork/tables/SLA_Supplementary_Tables.xlsx` | The workbook as approved and submitted, kept **only as the historical comparison baseline** that `make verify` checks against. Not the current supplementary material. |

The corrected workbook is built from the original, preserving its formatting,
sheet names, column order and number formats; every sheet is verified against its
regenerated CSV afterwards. Two supplementary table values were corrected —
Table S6's class I position 312 `coverage` and Table S4b's `class1_F_alt`
identity. Neither correction changes a figure or manuscript conclusion.
`CORRECTION_LOG.md` gives the evidence for each.

## Optional: recompute the structural measurements

The default workflow assembles Figures 5/6 from the archived rendered panels and
carries Tables S4a/S4b. To recompute those tables from the deposited model files:

```bash
micromamba env create -f environments/structures.yml
micromamba activate slahla-pub-structures
make structures && make tables && make verify
```

This needs PyMOL, because the published quantities are defined by PyMOL's
`cealign`, `super` and `align` — three different superpositions of the same pair,
plus groove and non-groove domain RMSDs. Substituting another aligner would not
reproduce them.

## Citing

See `CITATION.cff`. Cite the manuscript as the primary reference and this
repository by its release:

> Reproduction package for *Comparative Bioinformatic Analysis of Human and
> Swine MHC Sequences: Implications for Xenotransplantation*, release `v1.0.0`,
> https://github.com/olpmi/sla-hla-publication

**No archive DOI has been issued**; `PROVENANCE.md` records it as pending.

**Model checkpoints are not deposited.** The recovered seed-42 checkpoints exist
but have not been uploaded; seeds 43–46 are unavailable and cannot be reproduced
by retraining. The published findings rest on the archived predictions from all
five seeds, which are included here under `data/archived_predictions/`.

## Layout

```
src/slahla_pub/      the analysis and figure code
configs/             paths and settings; nothing is hard-coded to a machine
data/                archived scientific inputs (read-only)
publication_artwork/ the approved figures and tables, as the comparison baseline
validation/          the eplet mask selection rules and membership fingerprints
tests/               the retained test set
```

`LICENSE` (MIT) covers the code. It does **not** cover third-party data or model
weights; `DATA_LICENSES.md` states the terms for each.
