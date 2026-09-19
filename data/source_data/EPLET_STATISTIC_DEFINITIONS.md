# Three eplet quantities that are easy to conflate, and are not the same

`revised_eplet_statistics.csv` carries three numbers that all read in English as
"the eplet contrast". They weight the data differently and they answer different
questions. Which one is quoted changes the value, and the permutation test does
**not** test the number the headline quotes.

## 1. Covered-observation-weighted mean — the column `difference`

```
eplet_mean  = mean of the agreement indicator over EVERY (pair, position)
              observation at an eplet-associated position
other_mean  = the same over every observation elsewhere
difference  = eplet_mean - other_mean
```

Every observation counts once, so a position covered by 200 pairs contributes 200
times as much as a position covered by one. Because eplet-associated positions
are better covered than others (class I coverage 0.949 against 0.874; class II
0.756 against 0.580), this weighting is not neutral between the two groups.

Values, primary `agree_gap_skip` reading:

| class | eplet_mean | other_mean | **difference** | pair-bootstrap 95% CI |
|---|---|---|---|---|
| class I | 0.590 | 0.782 | **−0.192** | −0.196 to −0.187 |
| class II | 0.723 | 0.837 | **−0.113** | −0.122 to −0.104 |

## 2. Equal-position mean — the column `profile_difference`

```
for each HLA locus:
    collapse to one mean agreement per POSITION
    contrast = mean over eplet positions - mean over other positions
profile_difference = average of those per-locus contrasts, weighted by the
                     number of positions in the locus
```

Each position counts once within its locus regardless of how many pairs cover it,
so differential coverage cannot drive the contrast. Loci are weighted by how many
positions they contribute, not by how many pairs.

| class | **profile_difference** |
|---|---|
| class I | **−0.178** |
| class II | **−0.137** |

## 3. The permutation-test statistic — what `circular_shift_p` actually tests

The circular-shift test rotates the eplet mask within each locus and recomputes
the contrast. **Its observed statistic is `profile_difference`, the equal-position
mean — not `difference`.** The null distribution is built from the same
equal-position profile.

| class | statistic tested | value | p, reported as (extreme + 1)/(shifts + 1) |
|---|---|---|---|
| class I | `profile_difference` | −0.178 | < 1×10⁻⁴ |
| class II | `profile_difference` | −0.137 | 6.0×10⁻⁴ |

**Do not attach the permutation p-value to the observation-weighted
`difference`.** They are different statistics; the p-value belongs to the
equal-position contrast.

## What to quote

* Quote **`difference` with its pair-bootstrap interval** as the effect estimate,
  and say it is observation-weighted.
* Quote **`profile_difference` with `circular_shift_p`** as the test, and say the
  test statistic is the equal-position contrast.
* Never present one of the two means with the other's uncertainty or p-value.

The `mannwhitney_p_descriptive` column is a third thing again: a position-level
test whose independence assumption is violated here by repeated comparators,
within-pair dependence and spatial autocorrelation. It is **descriptive only** and
is not an inference.

## The DQB1 comparison stays descriptive

Class II eplet depletion is larger at DQB1 than at the other class II loci. **No
between-locus test has been performed**, so that comparison is reported as a
description of the per-locus profiles and not as a finding. A claim that DQB1
differs from the other loci would need a test that respects the same dependence
structure the primary analysis does, and none has been run.
