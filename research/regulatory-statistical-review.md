# Regulatory v1 statistical review

## Review status

The protocol and input/QC checks below are complete. Primary predictions, error
metrics, bootstrap draws and spatial-control contrasts reproduce independently.
The final exact-source comparator retains **26 failed secondary exact-rank
metrics**; a separately versioned post-hoc sensitivity study diagnoses the
remaining arithmetic fragility without changing their failed v1 status.

The executable frozen authority is
[regulatory-prediction-v1.json](../experiments/regulatory-prediction-v1.json),
SHA-256 `c8e97fabef2afe8ccaa5a01a7ef7bb952f9c9f15857e5c540af112a66e78cb58`.
The first draft was amended before new outcome acquisition to retain incomplete
positions and use chemistry-stratified rotations. The final input-selection and
model choices were not changed after acquisition.

## Independently checked eligibility

A separate direct CSV pass reproduced the frozen numeric/tag/reference criteria,
without fitting or inspecting positional outcome patterns:

| Primary element | Eligible SNVs | Represented positions |
| --- | ---: | ---: |
| F9 | 888 | 303 |
| FOXE1 | 1793 | 600 |
| GP1BB | 1154 | 385 |
| HBG1 | 820 | 274 |
| HNF4A | 854 | 285 |
| LDLR | 945 | 318 |
| MSMB | 1772 | 591 |
| PKLR | 1369 | 470 |
| TERT | 770 | 259 |
| **Primary total** | **10365** | **3485** |

Previously exposed HBB retains 561 SNVs at 187 positions. Validation-only records
retain LDLR repeat 950, PKLR later timepoint 1379, TERT other cell 770, TERT
knockdown 769 and TERT scrambled control 769.

All source positions remain represented, but not every alternate allele passes
QC. F9 has 17 two-allele and two one-allele positions; PKLR has 39 two-allele and
one one-allele positions. Dropping whole incomplete positions would discard
qualified data and change the stated target population. The final protocol keeps
eligible alleles and uses position-level error weighting.

Neither source p-values nor coefficient magnitudes determine inclusion. Source
deletion/equality rows do not enter these prediction fits. An original-source
crosswalk is required to interpret those representations separately.

## Guard and resampling checks

HBB physical test intervals are 1–37, 38–74, 75–112, 113–149 and 150–187.
Guarded training-position counts are respectively 130, 110, 109, 110 and 129.
Every fold's minimum training/test center distance is 21 bp.

An 11-bp motif's variant-overlapping windows can collectively reach ten bases in
either direction. Radius-ten supports are disjoint only when centers are more
than 20 bp apart. A ten-base center guard alone would leave overlapping sequence
support. The implemented contract therefore requires distance strictly greater
than 20.

The statistics helper tests establish deterministic SplitMix64 streams, common
resampling draws across the three paired comparisons, position weighting inside
partial physical blocks, Type-7 quantiles, Holm within-family arithmetic,
Bonferroni across all six planned tails, and plus-one tail counting with ties.
These arithmetic checks do not establish sampling independence or calibrated
inference.

Implementation conventions are explicit: one null stream over the union of HBB
and the nine primary loci, in full-URN then directed-mutation order; singleton
strata consume no draw. Each bootstrap family reuses one index table across all
three comparisons.

## Numerical artifacts found during review

Two implementation issues were identified before accepting reproduction:

1. A naive floating-point mean can make the computed variance of an exactly
   constant decimal vector slightly positive. Correlation must explicitly reject
   an exactly constant input before centering; a variance-greater-than-zero check
   alone is insufficient. A seven-element constant 0.1 vector is a regression
   control.
2. Even after that fix, averaging a constant prediction across groups containing
   one, two or three alleles can produce one-ULP differences between position
   means. Ranking then turns the differences into an artificial Spearman
   association. An initial independent HBB-to-F9 mean-baseline output showed
   undefined row correlation but position Spearman approximately −0.114, although
   every row prediction was the same constant. The reduction must preserve a
   mathematically constant group exactly before ranking.

Both owners were asked to retain the original attempts and regenerate corrected
artifacts. These are software/numerical corrections, not new features,
hyperparameters or biological findings. Primary squared errors are unaffected.

A subsequent full comparator found 130 remaining secondary Spearman differences:
the same ranking algorithm can assign different ranks when source-decimal means
or nearly tied model outputs differ by one ULP. Recomputing all 240 saved Rust
row/position correlations independently from **identical vectors** produced zero
failures (maximum difference about `5.27e-16`). The remaining issue is therefore
the numerical representation of ties, not a different correlation formula.

Exact source-decimal position means and reuse of identical-feature predictions
were implemented as semantic arithmetic corrections, reducing the mismatch count
from 130 to 26. The frozen tolerances were not widened and no rank epsilon was
added. The remaining v1 raw endpoint is explicitly failed.

The separately frozen
[numerical-rank sensitivity supplement](rank-sensitivity-review.md) evaluates
all five requested prediction-reporting grids plus raw exact ranks. All five
grids reproduce all 240 row/position metrics, whereas raw exact ranks reproduce
214. These are different, explicitly post-hoc rank estimands; no grid is selected
as a retroactive fix for v1 or evidence of genetic novelty.

## Primary loss and control results

These results do not depend on the unresolved secondary ranking issue. A third
artifact comparison matched **131112** prediction/observation/fold identities
across all four evaluation families and six models. Maximum prediction
difference was approximately `7.73e-12`, below the frozen `1e-7` tolerance.
Both 2000-draw bootstrap index tables matched exactly; maximum paired bootstrap
statistic difference was approximately `9.03e-15`. All six sets of 99 rotation
contrasts matched to about `1.38e-14`.

### Fixed-model errors

Position-weighted MSE, in squared log2-effect units:

| Model | Guarded HBB | Primary nine-element holdout | HBB row-random diagnostic |
| --- | ---: | ---: | ---: |
| mean | 0.021934 | 0.116392 | 0.021462 |
| chemistry | 0.022683 | 0.116899 | 0.021726 |
| coordinate | 0.026120 | 0.116529 | 0.020792 |
| delta3mer | 0.025027 | 0.118371 | 0.020942 |
| motif | 0.024997 | 0.116829 | 0.021855 |
| full | 0.029610 | 0.117989 | 0.020436 |

The primary result is negative: no fitted feature model improves on the
training-mean baseline in either guarded HBB or the nine-element macro
evaluation. The full model's favorable row-random diagnostic does not survive
the guarded split. This is a split/complexity warning, not a causal estimate of
leakage or a biological mechanism.

### Every prespecified comparison

Improvement is baseline MSE minus candidate MSE; positive is favorable. Intervals
are descriptive conditional resampling ranges, not population guarantees.

| Family | Comparison | Improvement | Bootstrap 95% range | Rotation tail | Holm within three | Bonferroni across six |
| --- | --- | ---: | --- | ---: | ---: | ---: |
| Guarded HBB | delta3mer vs chemistry | −0.0023442 | [−0.0043825, +0.0001200] | 0.08 | 0.24 | 0.48 |
| Guarded HBB | motif vs chemistry | −0.0023135 | [−0.0039156, −0.0007471] | 1.00 | 1.00 | 1.00 |
| Guarded HBB | full vs coordinate | −0.0034899 | [−0.0056994, −0.0007149] | 0.27 | 0.54 | 1.00 |
| Primary holdout | delta3mer vs chemistry | −0.0014717 | [−0.0025342, −0.0001316] | 0.98 | 1.00 | 1.00 |
| Primary holdout | motif vs chemistry | +0.0000697 | [−0.0005804, +0.0006429] | 0.01 | 0.03 | 0.06 |
| Primary holdout | full vs coordinate | −0.0014601 | [−0.0025373, −0.0001647] | 0.95 | 1.00 | 1.00 |

The small motif-versus-chemistry difference is about 0.06% of chemistry MSE,
has an element-resampling interval spanning zero, and does not beat the simpler
mean model. Selecting the favorable within-family tail while omitting the other
family or simpler baseline would overstate the result. The artificial rotation
tail remains a robustness diagnostic even after multiplicity arithmetic.

### Context-blind finite bounds

For two observed coefficients `y1,y2` and one context-blind prediction `a`:

`((y1-a)^2+(y2-a)^2)/2 = (a-(y1+y2)/2)^2+(y1-y2)^2/4`.

Thus the observed paired MSE has lower bound one quarter of the mean squared
coefficient difference. This algebra does not identify its biological cause.

| Source-defined pair | Matched QC-passing SNVs | Paired coefficient MSE | Context-blind lower bound |
| --- | ---: | ---: | ---: |
| LDLR library repeat | 943 | 0.0339334 | 0.0084834 |
| PKLR 48h versus 24h | 1368 | 0.0569433 | 0.0142358 |
| TERT SF7996 versus HEK293T | 770 | 0.1028804 | 0.0257201 |
| TERT GABPA knockdown versus scrambled control | 769 | 0.0338541 | 0.0084635 |
| TERT scrambled versus unperturbed SF7996 | 769 | 0.0277620 | 0.0069405 |

Each row is its own paired dataset. Reused TERT observations are not pooled as
independent repeats. The LDLR comparison shows that library/measurement
differences alone can create a nonzero bound. Different genes, normalization,
libraries, interventions and cells prevent interpreting a larger bound as a
clean causal estimate of cell-state variation.

## Interpretation requirements for the final artifacts

- The primary outcome remains position-weighted, equal-element macro MSE.
  Secondary rank correlations must not replace it after seeing unfavorable
  results.
- Report all three prespecified candidate/baseline comparisons, all nine
  held-out elements, both control families and all model attempts. A favorable
  single element or family is not the primary result.
- Mutation-stratified rotations preserve class-specific score distributions and
  spatial ordering, but break some cross-class correlation and do not establish
  promoter stationarity. Their empirical tails and multiplicity adjustments
  remain controls, not calibrated biological-significance p-values.
- Position and element bootstrap summaries condition on frozen out-of-fold
  predictions. They omit fitting/selection uncertainty and cannot turn nine
  selected, same-study elements into independent laboratories.
- Guarded and row-random folds differ in both leakage opportunity and split
  difficulty. Their performance gap is a diagnostic, not a causal estimate of
  leakage.
- Context-blind variance bounds must use unique matched variant/context
  observations within each stated comparison. Reused TERT contexts must not be
  counted repeatedly as independent observations when aggregating pairs.
- Source normalization, library construction, cell state and assay noise can all
  contribute to context differences. A finite lower bound is algebra, not a
  causal allocation of those differences.
- Agreement of independent Rust and Python predictions checks numerical
  reproduction of the shared specification. It does not prove the specification
  describes biology or that both avoid a shared source-assumption error.
