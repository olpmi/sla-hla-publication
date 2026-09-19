# Reproducibility

What the published results start from, how to rerun the stages upstream of that,
and what is known not to reproduce.

## Two validations, and the difference between them

This repository has been run two ways, and the reports distinguish them:

- **With locally available inputs** — the eplet position mask supplied from a
  registry export held under licence. All 16 figures and all publication tables
  are produced, and the eplet statistics are recomputed.
- **From the public bundle alone** — no mask. 13 of 16 figures; Figures S8, S9
  and S10 are not produced, and Figure 4 and Table S3 fall back to archived
  aggregate statistics. The run reports this explicitly and records it in
  `outputs/reports/missing_inputs.json`.

`VALIDATION_REPORT.md` states which mode each result came from. Nothing in the
public-bundle run is skipped silently, and no output is produced from a guessed
substitute for a missing input.

## Scientific starting points

`make reproduce` does not re-derive everything from raw sequence. It starts from
archived scientific inputs and recomputes forward. Each starting point is named
here and in `PUBLICATION_OUTPUTS.csv`.

| Stage | Starts from | Recomputes |
|---|---|---|
| Classifier | the **ten archived five-seed prediction tables** (`data/archived_predictions/`, 5 seeds × 2 classes) and their run configs | ensemble probabilities, per-sequence stability, cross-seed agreement, leave-one-seed-out ensembles, held-out validation summary |
| Eplet positions | `data/eplet_agreement/identity_eplet_agreement.parquet` — 112,820 (pair, position) agreement observations under both gap readings — **plus the eplet position mask, which is not shipped and must be supplied locally** | observation-weighted and equal-position contrasts, pair-bootstrap intervals, within-locus circular-shift permutation |
| Sequence / haplotypes | the curated haplotype workbooks in `data/haplotypes/` | the haplotype half of Table S1 |
| Structural (default) | the **archived rendered PyMOL panels** in `data/structural_panels/` | Figures 5/6 panel assembly, label projection and placement |
| Structural (optional) | the **24 deposited model PDBs** and 3 experimental reference structures | Tables S4a and S4b in full |

Where a stage starts from an archived upstream result rather than raw input, it
says so. The rendered panels in particular are an archived result: the default
workflow re-does the cropping, projection, label placement and annotation, but
does not re-render the rasters.

## Provenance of the classifier archive

The published classifier findings rest on the **archived original predictions**.
These are primary and are what ships here.

Checkpoint availability is a separate matter from prediction availability:

- **Seed 42, both classes** — the original weights were recovered. Reloading them
  reproduces the archived probabilities to a maximum absolute difference of
  1.34 × 10⁻⁵ (class I) and 1.82 × 10⁻⁵ (class II); that is inference
  nondeterminism, not training variation.
- **Seeds 43–46, both classes** — the original weights are **unavailable** and
  cannot be reproduced. Training did not force deterministic GPU kernels, so
  retraining under the identical protocol yields different weights. Anyone who
  retrains is producing a *new* draw from the procedure, not recovering the
  archived one.
- Reconstructed checkpoints produced by later retraining exist in the private
  development record and are **not** interchangeable with the originals. They are
  excluded here. No reconstructed prediction is substituted for an archived one
  anywhere in this repository, and `slahla_pub.ensemble` refuses to load anything
  outside the canonical ten.

The archive also contains two controlled split-design arms (`armA_full_random`,
`armB_dedup_random`, seed 42) built for the reviewer response, and two legacy
single-model tables from the submitted analysis. None is an ensemble member;
`ensemble.probs_path` raises if asked for one.

## Rerunning the upstream stages

None of this is needed to reproduce the published results.

### Inputs and how to acquire them

`data/MANIFEST.json` records 54 third-party files with source URL, SHA-256, byte
count, retrieval timestamp and release.

| Resource | Release | Acquisition |
|---|---|---|
| IPD-IMGT/HLA protein sequences and alignments | **v3.65.0-alpha** | `raw.githubusercontent.com/ANHIG/IMGTHLA/v3.65.0-alpha/...` as recorded in the manifest |
| IPD-MHC (SLA) protein sequences | pinned by SHA-256 only | `ftp.ebi.ac.uk/pub/databases/ipd/mhc/MHC_prot.fasta` |
| CIWD 3.0.0 | 2020-03-20 | 18th IHIW distribution, per the manifest URLs |
| HLA Eplet Registry export | accessed 2026-08-05 | obtain directly from https://www.epregistry.com.br/ (account required) |
| Experimental reference structures | — | 1AQD, 1HHK, 1S9V; **included** here (RCSB PDB archive data are CC0) |

Verify every download against the manifest checksum before use. See
`DATA_LICENSES.md` for what may and may not be redistributed.

### Software

The published analyses ran under Python 3.11.15 with PyTorch 2.13.0,
Transformers 4.57.6, scikit-learn 1.9.0, NumPy 2.4.6 and Biopython 1.87.
External programs: DIAMOND ≥ 2.1.9 (all-vs-all protein alignment, e-value
threshold 1 × 10⁻³), MAFFT ≥ 7.5 (`--keeplength --add`, projecting SLA onto HLA
registry coordinates), PyMOL ≥ 3.0 (superposition and rendering), ColabFold over
AlphaFold2 (structure prediction). `environments/upstream.yml` installs the
CPU-side tools; PyTorch with CUDA must be installed to match your GPU.

### Classifier training

Settings are in `configs/upstream.yaml`: `Rostlab/prot_bert_bfd`, max length 512,
5 epochs, batch size 8, learning rate 1 × 10⁻⁵, weight decay 0.01, StepLR
(step 3, γ 0.1), one masked residue at each terminus, exact-sequence
deduplication, 99 %-identity clusters held out whole, validation fraction 0.2,
seeds 42–46.

**Checkpoint selection is the final epoch.** There is no early stopping and no
best-epoch selection; the run config each member ships records the final
validation accuracy and loss.

Compute: the canonical queue of 14 runs took about 9 hours on an A100. Training
on CPU is not practical.

### Expected differences from the archive

A rerun is a new draw. Acceptance criteria were fixed in advance of any
retraining and must not be widened afterwards:

| | Criterion |
|---|---|
| Split identity | **exact**. The split is a seeded deterministic function of the record set, not of the training run. A mismatch invalidates the comparison before any weight is examined. |
| Validation accuracy | within 0.0336 (class I) / 0.0122 (class II) of the recorded value — twice the observed across-seed SD of the archived runs. |
| Per-sequence locus agreement | ≥ 0.90 against the archived same-seed table. Subgroup agreement is deliberately *not* a criterion: across the five archived seeds no class I SLA sequence keeps the same subgroup label, so requiring it would fail a correct rerun. |
| Probability agreement | reported, never thresholded. Retraining on different hardware will not reproduce archived probabilities to 1 × 10⁻⁶ and nothing here assumes it will. |
| Ensemble stability | ≥ 0.95 agreement on the ensemble top locus; ≤ 0.05 per-locus mass difference. The ensemble is the reported estimator, so it is the ensemble that must be stable, not any member. |

## Known limits of the published results

These affect the published findings and are not development clutter.

**Gap handling is a stated choice.** Under per-sequence pairing many SLA records
are partial, and eplet-associated positions are systematically better covered
than others (class I 0.949 vs 0.874; class II 0.756 vs 0.580). Scoring gaps as
mismatches confounds coverage with divergence and, for class II DRB1, **reverses
the sign** of the contrast. Gapped positions are therefore excluded in the
primary reading, and both readings are reported. Every contrast quoted anywhere
must name its gap rule.

**Two eplet statistics are easy to conflate.** `difference` is
observation-weighted: every (pair, position) observation counts once, so a
position covered by 200 pairs counts 200 times. `profile_difference` is
equal-position: one mean per position within each locus, averaged across loci by
position count. **The circular-shift p-value tests `profile_difference`, not
`difference`.** Never present one with the other's uncertainty.

**Position-level Mann–Whitney and t-tests are descriptive only.** Their
independence assumption is violated here: one HLA comparator serves up to 151 SLA
sequences, and hundreds of positions come from one alignment. Primary inference
is the pair-level bootstrap plus the within-locus circular-shift permutation.

**No between-locus test was performed.** Class II eplet depletion is larger at
DQB1 (+14.9 percentage points excess mismatch) than DRB1 (+9.0) under the primary
reading, but that comparison is a description of the per-locus profiles, not a
finding. A claim that they differ would need a test respecting the same
dependence structure, and none has been run.

**Subgroup-level assignment is withdrawn.** Across the five canonical models, **no
class I SLA sequence (0 of 260)** and only **4 of 152 class II sequences** keep
the same HLA subgroup label. Ensembling does not recover it. Only locus-level
probability mass is reported, with its leave-one-seed-out range. Note that the
neighbouring quantity `same_top_subgroup_all_loo` (0.1654 class I, 0.5263 class
II) measures agreement among the five *leave-one-seed-out ensembles* and is a
different endpoint; quoting it as model unanimity overstates stability roughly
sixteenfold.

**Class I accuracy is over 53 of 54 labels.** `C_Other` (26 sequences in a single
99 %-identity cluster) has no validation side under the cluster split. The label
count travels with every accuracy figure.

**DIAMOND reports only alignments passing its threshold**, so low-similarity
comparisons may be absent and reported locus means can be biased upward.

**The structural models are monomers.** They come from ColabFold over AlphaFold2,
not the AlphaFold 3 Server named in the submitted Methods. A ColabFold monomer
job is one chain, so the class I models carry no β2-microglobulin and no peptide,
and the class II models no α chain and no peptide.

**The final figures' RMSDs are reproducible; the submitted ones are not.** Each
published panel prints `round(rmsd_cealign_CA, 1)` of the whole-chain Cα CE RMSD,
and all nine agree with the recomputed measurement, the render report and the
final structural manifest — see `outputs/reports/structural_panel_audit.csv`.
Figure 5F prints 2.6 Å for HLA-C\*07:02, which is the recomputed 2.638 Å rounded;
there is no discrepancy there.

What has no computational provenance is the *submitted* version's printed values,
which were produced interactively and typed in. They are kept separately in
`TableS4b_provenance_submitted_version.csv` — including the submitted legend's
HLA-C\*07:01 and its 3.4 Å — so they cannot be read as current values. Table S4b
retains the three named superposition methods (`cealign`, `super`, `align`) with
their units and denominators, because they answer which quantity a given number
is. The `class1_F_alt` row superposes the C\*07:01 model the submitted legend
named; it is labelled a historical comparison and is not a published panel.

**The cohort is a dated snapshot.** 452 retained SLA protein records, of which 260
class I and 152 class II enter the classifier analysis. Thirty-one alleles have
been withdrawn or renumbered since the submitted analysis, whose input files no
longer exist. The revised analysis is pinned to a checksummed distribution and
does not claim to use the same records.
