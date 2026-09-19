# Correction log

Every value the regenerated outputs change relative to the approved publication
materials, with the evidence and the downstream effect. The approved copies under
`publication_artwork/` are **not** modified — they are the comparison baseline
`make verify` checks against — so the corrected values appear only under
`outputs/`.

`make verify` reports C1 and C2 as corrections and fails on any difference not on
this list.

## C1 — Table S6, class I position 312, `coverage`

| | |
|---|---|
| Approved | `0.0454545454545454` |
| Corrected | `0.003861003861003861` |
| Calculation | `n_pairs_covering / class pair count` = 1 / 259 |
| Source | `data/eplet_agreement/identity_eplet_agreement.parquet` |
| Fixed in | `slahla_pub.tables.table_s6` |

The approved transformation divided `n_pairs_covering` by the number of
agreement-matrix rows present at that position. For every position covered by the
whole cohort that equals the class pair count, so 526 of 527 rows were right by
coincidence. Class I position 312 appears only in the 22 HLA-A pair rows, so it
was divided by 22. The denominator is now a fixed per-class constant
(`S6_PAIR_DENOMINATOR` = 259 class I, 151 class II), asserted against the pair
count in the matrix before use, so the error cannot reappear at another position.

**Downstream effect: none.** Table S6 has no rank, priority or score column and is
ordered by position, so nothing reorders. `coverage` appears in no figure, no
reported statistic and no scientific endpoint check. Figures S8 and S9 carry a
*different* quantity, `coverage_share_of_pairs` — coverage within each comparator
locus, 1/22 at this position — which is correct as published and unchanged.

## C2 — Table S4b, `class1_F_alt`, `percent_identity_over_shared_positions`

| | |
|---|---|
| Approved | `74.92` |
| Corrected | `75.58` |
| Calculation | 229 identical residues / 303 positions where both chains have a residue |
| Method | MAFFT v7.526 (2024/Apr/26), `mafft --quiet --auto` |
| Evidence | `data/carried/tableS4b_pair_identity_verified.csv` |
| Fixed in | `slahla_pub.structural.pair_identity` |

The approved table repeats the `class1_F` value on the `class1_F_alt` row: 74.92%
is the identity of the C\*07:02 model against SLA-3\*06:01, not of the C\*07:01
model. Recomputed independently from the two modelled sequences:

| | `class1_F` | `class1_F_alt` |
|---|---|---|
| HLA model | `HLA_C07_C_07_02_unrelaxed_rank_001_…pdb` | `C07_C0701_unrelaxed_CUSTOM_B.pdb` |
| HLA model SHA-256 | `3b82887c61b0f463…` | `7334604f12f16147…` |
| SLA model | `SLA_C07_SLA-3_06_01_unrelaxed_rank_001_…pdb` | same |
| SLA model SHA-256 | `22e199ccf66db8ef…` | `22e199ccf66db8ef…` |
| modelled residues | 307 / 303 | 309 / 303 |
| identical residues | 227 | 229 |
| positions both present | 303 | 303 |
| identity | 227/303 = **74.92%** | 229/303 = **75.58%** |

The nine other rows reproduce the approved value exactly under the same
procedure, which is what makes the tenth a defect rather than a method
difference.

**Downstream effect: none.** `class1_F_alt` is a historical comparison against the
allele the submitted legend named. It is not a displayed panel, and `record_level`
says so. **Figure 5F is unchanged and correct**: it shows HLA-C\*07:02 and prints
2.6 Å, which is `round(2.638, 1)` of its own recomputed CE RMSD.

Without MAFFT installed both identity columns fall back to the approved values,
and `identity_method` records that.

## Table S1 — two notes on the current output

**`local_name` agrees on all 82 rows.** An earlier parser in this repository took
`local_name` from the first non-empty value across every haplotype an allele
appeared in, which imported a miniature-line label (X, Y, Z) onto a commercial
haplotype. `slahla_pub.haplotypes.alleles_by_haplotype` now follows the selected
haplotype row. The approved table was correct and needed no change on this
account.

**Three haplotype selections are carried, not recomputed.** `SLA-DQA*01:01`,
`SLA-DQA*02:02:02` and `SLA-DRA*02:03:01` each sit in two curated haplotypes and
the approved table reports one. Four candidate rules — first row, last row, first
commercial row, last commercial row — were each tested against all 34
multi-haplotype alleles and each fails on at least one, so the rule is not
recoverable from the shipped inputs. The archived selection is carried, recorded
in `data/carried/tableS1_haplotype_selection.csv` with its source, and labelled in
the `haplotype_selection_source` column. The full regenerated set is retained
alongside in `haplotypes_all_curated`. No manuscript result depends on these
accessory-locus rows.
