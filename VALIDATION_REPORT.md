# Validation report

Run on 2026-09-19 · Python 3.11.15 · macOS-26.6.2-arm64-arm-64bit

numpy 2.4.6 · pandas 3.0.5 · scipy 1.17.1 · matplotlib 3.11.1 · pillow 12.3.0 · PyMOL present

## What was executed

```
make reproduce      # ensemble -> eplet statistics -> tables -> figures -> verify
make structures     # optional, needs PyMOL: recompute Tables S4a and S4b
make test
```

## What matched

**Publication tables** — 7 of 7 identical to the approved materials on every column and row: S1, S2, S3, S4a, S4b, S5, S6.

**Recomputed intermediates** — 7 of 7 identical to the archived tables they were recomputed from:

- `Figure4_plotted_values.csv`
- `canonical_ensemble_per_sequence_stability.csv`
- `canonical_ensemble_predictions.csv`
- `canonical_seed_stability.csv`
- `canonical_validation_across_seeds.csv`
- `eplet_statistics_v6.csv`
- `leave_one_seed_out_by_locus.csv`

**Scientific endpoints** — 47 of 47 checks passed. These are the manuscript's own endpoints with their definitions attached, including the cross-seed subgroup unanimity (0/260 class I, 4/152 class II), the gap-rule-labelled eplet contrasts, and the distinction between the observation-weighted `difference` and the equal-position `profile_difference` that the permutation test actually tests.

**Stochastic quantities reproduced exactly.** The pair-bootstrap intervals (2000 resamples) and circular-shift p-values (10,000 shifts) match the published values digit for digit, because the original seed (20260913), resample counts and sampling procedure were recovered and re-executed. No tolerance was applied to them.

**Inputs unchanged.** `data/source_data/` (85 files) is byte-identical after a full run: confirmed.

## Figures

Compared on content and layout, not bytes: PDFs embed creation timestamps and font subsets, so byte equality would fail on identical artwork.

| figure | canvas matches | mean abs. pixel difference | ink bbox shift | PDF | SVG |
|---|---|---|---|---|---|
| Figure1 | no | 0.04225 | 4 px | yes | yes |
| Figure2 | no | 0.03138 | 7 px | yes | yes |
| Figure3 | no | 0.04606 | 4 px | yes | yes |
| Figure4 | no | 0.01615 | 1 px | yes | yes |
| Figure5 | no | 0.01488 | 2 px | yes | yes |
| Figure6 | no | 0.00846 | 1 px | yes | yes |
| FigureS1 | yes | 0.00845 | 3 px | yes | yes |
| FigureS2 | no | 0.0269 | 2 px | yes | yes |
| FigureS3 | yes | 0.00336 | 2 px | yes | yes |
| FigureS4 | yes | 0.00952 | 2 px | yes | yes |
| FigureS5 | yes | 0.00494 | 1 px | yes | yes |
| FigureS6 | no | 0.02835 | 2 px | yes | yes |
| FigureS7 | yes | 0.01636 | 2 px | yes | yes |
| FigureS8 | no | 0.02348 | 2 px | yes | yes |
| FigureS9 | no | 0.02292 | 2 px | yes | yes |
| FigureS10 | no | 0.02313 | 2 px | yes | yes |

All sixteen figures regenerate in PNG, PDF and SVG. The residual pixel difference is text antialiasing: the approved rasters were rendered by a different freetype build, so glyph edges differ by a fraction of a pixel while the plotted geometry does not. Multi-panel figures are cropped to their ink (`bbox_inches='tight'`), which makes the canvas a few pixels wider or narrower when glyph metrics shift; the five single-panel supplementary figures use the full canvas and match it exactly.

## Structural panels: measured, displayed, reconciled

All 9 displayed panels agree. The value each panel prints is `round(rmsd_cealign_CA, 1)` of the whole-chain C-alpha CE RMSD, and three independent records of that measurement - the render report the figure builder reads, the final structural manifest, and the recomputation from the deposited model files - agree to within 0.0005 A.

| figure | panel | HLA | SLA | recomputed CE RMSD | displayed | agree |
|---|---|---|---|---|---|---|
| Figure 5 | A | A*02:01 | SLA-1*02:01 | 2.681 A | 2.7 A | yes |
| Figure 5 | B | A*03:01 | SLA-2*01:02 | 3.782 A | 3.8 A | yes |
| Figure 5 | C | B*08:01 | SLA-3*01:01 | 2.738 A | 2.7 A | yes |
| Figure 5 | D | B*15:01 | SLA-1*01:01 | 2.966 A | 3.0 A | yes |
| Figure 5 | E | B*57:01 | SLA-2*01:01 | 2.413 A | 2.4 A | yes |
| Figure 5 | F | C*07:02 | SLA-3*06:01 | 2.638 A | 2.6 A | yes |
| Figure 6 | A | DRB1*01:02 | SLA-DRB1*01:01 | 3.102 A | 3.1 A | yes |
| Figure 6 | B | DRB1*11:04 | SLA-DRB1*04:02 | 3.374 A | 3.4 A | yes |
| Figure 6 | C | DQB1*02:02 | SLA-DQB1*02:02 | 1.44 A | 1.4 A | yes |

**Figure 5F carries no discrepancy.** The final legend and the final artwork both name HLA-C*07:02, and the panel prints 2.6 A, which is the recomputed 2.638 A rounded. The HLA-C*07:01 allele and the 3.4 A value belong to the *submitted* version of the manuscript. They are retained only as history, in `TableS4b_provenance_submitted_version.csv`, and the `class1_F_alt` row that superposes the C*07:01 model is labelled a historical comparison, not a displayed panel.

## Recomputed, carried, or not executed

Exact agreement with an archived table is not the same claim as independent recomputation. Every column of every publication table is classified and the classification is checked, not asserted: each `recomputed` column is recalculated here and compared. Full detail in `outputs/reports/column_provenance.csv`.

| class | columns | meaning |
|---|---|---|
| recomputed | 175 | recalculated in this run from packaged inputs |
| carried | 19 | archived value copied into an assembled output |
| external_not_executed | 12 | needs an input this repository does not ship |
| key | 11 | identifier column |

| table | recomputed | carried | external | key |
|---|---|---|---|---|
| S1 | 6 | 1 | 7 | 1 |
| S2 | 28 | 0 | 0 | 0 |
| S3 | 26 | 0 | 0 | 0 |
| S4a | 35 | 0 | 0 | 0 |
| S4b | 37 | 0 | 0 | 7 |
| S4b_provenance | 0 | 12 | 0 | 0 |
| S5 | 38 | 0 | 0 | 0 |
| S6 | 5 | 6 | 5 | 3 |

## Differences from the approved tables, and what causes them

**Table S1 - 0 differing cells of 492 compared.** The seven haplotype columns are regenerated from the curated workbooks in `data/haplotypes/`; all 82 alleles are covered, and `sla_locus`, `mhc_class`, `herd` and `local_name` match on every row. Cell-level detail is in `outputs/reports/tableS1_cell_differences.csv`.

All remaining differences are the same kind: an allele that appears in more than one curated haplotype, where the approved table reports one and this regeneration reports the first (`haplotype`) or all of them (`haplotypes`). The rule the approved table used to choose is not recorded and is not recoverable from the shipped inputs. The underlying workbook rows are identical either way, all three alleles are accessory-locus records (SLA-DQA, SLA-DRA) that enter no classifier, eplet or structural analysis, and no manuscript result depends on them.

A further five differences reported in the previous release have been **resolved**: `local_name` was being taken from the first non-empty value across every haplotype an allele appeared in, which imported a miniature-line label (X, Y, Z) onto a commercial haplotype. It now follows the selected haplotype row, and all 82 rows match.

**Table S1's identity columns** - `n_ciwd_hits`, `nearest_*`, `most_distant_*` - are pairwise identities against the CIWD 3.0 common-and-well-documented allele set. They need the DIAMOND all-vs-all stage over IPD-IMGT/HLA plus the CIWD release. That calculation was **not executed** in the clean-environment run and its inputs are not shipped; the values are carried from the approved table.

**Table S4b's two identity columns** are recomputed, not carried. The data dictionary states their method - MAFFT identity between the two modelled sequences, over positions where both chains have a residue - and nine of the ten rows reproduce exactly under it. The tenth, `class1_F_alt`, does not: the approved table repeats `class1_F`'s value (74.92% over 303 positions) for the alternative C*07:01 model, where MAFFT gives 75.58%. `class1_F_alt` is a historical comparison row, not a published panel, and no manuscript claim uses it. The approved value is left in place and the difference recorded rather than silently corrected. Without MAFFT both columns fall back to the approved values and `identity_method` says so.

**Table S6** is a mixed table. Its agreement half - `mean_agreement`, `divergence`, `n_pairs_covering`, `is_eplet_associated` - is recomputed here from `identity_eplet_agreement.parquet` and the eplet position mask, and matches the approved table exactly (max difference 1.1e-16). `coverage` matches for 526 of 527 rows; class I position 312, which a single pair covers, uses a 22-pair denominator in the approved table where every other row uses 259. That is recorded, not altered. The eplet-name half - `eplet_names`, `eplet_loci_defining`, `eplet_name_source`, `eplet_named_in_registry_group_mask`, `eplet_analysis_support` - requires the HLA Eplet Registry export, which cannot be redistributed; to rebuild it, obtain the export from https://www.epregistry.com.br/, pool the flattened `resi` column per locus, and join on `position`. The structural columns are carried from `structural_candidate_annotations_v6.csv`.

**The original seed 43-46 checkpoints** are unavailable, and retraining cannot reproduce them: GPU kernels were not made deterministic. The archived predictions those checkpoints produced are what the published classifier findings rest on, and those are shipped and recomputed from.

**The RMSDs printed in the SUBMITTED figures** were produced interactively and typed in; no computational provenance survives for them. They are retained in `TableS4b_provenance_submitted_version.csv` as history. The values the FINAL figures print are reproducible and reproduced: see the panel audit above.

## Release validation: two modes

The candidate was run two ways, because what it can produce depends on an input it does not ship. Both are reported; neither is presented as the other.

### A - with locally available inputs

the candidate with every input available locally: the eplet position mask generated by the repaired importer from the pinned registry export, plus PyMOL and MAFFT.

| command | exit |
|---|---|
| `make structures` | 0 |
| `make reproduce` | 0 |
| `make test` | 0 |

Figures produced: **16 of 16**. Column provenance: 175 recomputed, 19 carried, 12 external_not_executed, 11 key. Scientific checks: 47/47 passed. Tests: 113 passed.

Corrected supplementary workbook vs the regenerated CSVs: **AGREE**.

### B - from the public bundle alone

the public bundle alone, extracted to a fresh directory: no eplet position mask, no PyMOL, no MAFFT, no access to the private repositories.

| command | exit |
|---|---|
| `pip install -e . --no-deps` | 0 |
| `make reproduce` | 0 |
| `make test` | 0 |

Figures produced: **13 of 16**. Column provenance: 105 recomputed, 88 carried, 13 external_not_executed, 11 key. Scientific checks: 47/47 passed. Tests: 112 passed, 1 skipped (the skipped test requires the mask).

Not produced: **FigureS8, FigureS9, FigureS10** — each marks eplet-associated positions individually and the mask is absent. Produced from archived aggregate statistics instead of recomputed: Figure4, TableS3. the HLA eplet position mask, named in the MISSING EXTERNAL INPUTS block and recorded in outputs/reports/missing_inputs.json.

In both modes `data/source_data` is byte-unchanged (true) and the clean extraction has none in either mode's clean extraction. Optional environments: `environments/structures.yml` solves; pymol-open-source 3.1.0, mafft 7.526; `environments/upstream.yml` solves; diamond 2.2.7, mafft 7.526.

## Corrections applied in this release

| # | where | old | new | basis | downstream effect |
|---|---|---|---|---|---|
| C1 | Table S6, class I position 312, `coverage` | 0.0454545454545454 | 0.003861003861003861 | 1/259; the denominator is now a fixed per-class constant, asserted against the matrix | none - no rank column, no figure, no statistic |
| C2 | Table S4b, `class1_F_alt`, identity | 74.92% | 75.58% | MAFFT v7.526, 229/303 identical residues; evidence in `data/carried/tableS4b_pair_identity_verified.csv` | none - a historical comparison row, not a displayed panel. Figure 5F unchanged |
| C3 | Table S1, `local_name` (previous release) | miniature-line labels | the selected haplotype's value | a defect in this repository, not in the approved table | none - all 82 rows now agree |

Full detail, including the model checksums behind C2, is in `CORRECTION_LOG.md`. The approved copies under `publication_artwork/` keep their original values so the comparison baseline stays what the authors approved; `make verify` reports C1 and C2 as corrections and fails on anything not on that list.

## Table provenance

- **S1** — partial: haplotype columns regenerated from the curated workbooks; n_ciwd_hits/nearest_*/most_distant_* archived (need the upstream DIAMOND identity stage over IPD-IMGT/HLA and CIWD 3.0)
- **S2** — recomputed
- **S3** — recomputed
- **S4a** — recomputed
- **S4b** — recomputed except percent_identity_over_shared_positions and n_positions_both_present, which are archived
- **S5** — recomputed
- **S6** — archived (needs the licensed HLA Eplet Registry)

