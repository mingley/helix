# Rank resolution: a post-hoc numerical methods result

## Finding and scope

The frozen v1 exact-rank endpoint **remains failed**: 26 of 240 paired secondary
Spearman metrics differ beyond `1e-8`. The implementations nevertheless reproduce
131112 predictions within about `7.73e-12`, and agree on all primary errors,
folds, fits, bootstrap draws and spatial-control statistics.

This supplement finds that every prespecified finite prediction grid reproduces
all 240 correlations exactly. It diagnoses a tie-resolution problem, not a new
biological association. No model was refit, no p-value calculated, and no quantum
selected.

The separate pre-diagnostic
[protocol](../experiments/numerical-rank-sensitivity-v1.json) has SHA-256
`b968f0b3519b3890c36d2756ce97780f3a42effb8ee06b138e4cf395138b89ed`.
The parent analytical protocol remains unchanged at
`c8e97fabef2afe8ccaa5a01a7ef7bb952f9c9f15857e5c540af112a66e78cb58`.

## Every declared resolution

Predictions only are quantized after each implementation's native reduction.
Observed position means remain exact rational means of original source-decimal
scores. Integer grid codes use exact saved-float arithmetic and nearest
ties-to-even rounding.

| Reporting level, log2-effect units | Matched metric pairs | Mismatched | Both undefined | Maximum cross-implementation rho difference |
| --- | ---: | ---: | ---: | ---: |
| Raw exact | 214 / 240 | 26 | 36 | 0.007978226 |
| `1e-12` | 240 / 240 | 0 | 36 | 0 |
| `1e-10` | 240 / 240 | 0 | 36 | 0 |
| `1e-8` | 240 / 240 | 0 | 36 | 0 |
| `1e-6` | 240 / 240 | 0 | 36 | 0 |
| `0.01` | 240 / 240 | 0 | 36 | 0 |

There are zero undefined-status mismatches at every level. The 36 jointly
undefined cases remain separately visible rather than counted as defined
correlation estimates.

Agreement does not mean the coarser estimands have the same values as raw ranks.
Across the four finer grids, maximum change from either implementation's own
raw correlation is about 0.00417237. At 0.01, that maximum becomes about 0.108935.
The coarse level discards prediction-order detail. It would be misleading to
pick it because the implementations agree.

## Concrete worst case

For chemistry-only HBB-to-F9 transfer at the SNV-row level:

- Rust raw Spearman: 0.0619531934.
- Independent Python raw Spearman: 0.0539749674.
- Absolute difference: 0.0079782260.
- Maximum difference between their underlying predictions:
  `2.67e-16` log2-effect units.

The raw models each expose 12 numeric prediction values for the 12 mutation
identities. Mathematically tied fitted class effects can be separated differently
by final arithmetic operations. Exact ranking treats that separation as an
ordering, regardless of its negligible magnitude. A small absolute prediction
error therefore does not guarantee closeness of a discontinuous rank statistic
at ties.

One exact source symmetry makes this concrete: the HBB A>G and A>T training
classes each contain 36 rows and each has decimal-score sum 2/5, hence mean 1/90.
Their identical counts, means and fixed ridge treatment imply tied chemistry
class effects mathematically; tiny solver/reduction differences can split them
numerically. This is an arithmetic property of this selected source table, not
an assertion that those mutation classes are biologically equivalent.

The same-input check independently recomputed all 240 Rust correlations from
identical saved vectors with zero failures and maximum difference about
`5.27e-16`. The disagreement is not a different Spearman formula.

## Source precision and arithmetic corrections

The source CSV snapshot contains 16546 finite coefficients printed with two
fractional digits and 1922 with one. These are formatting facts about all
retained raw tables, not a measurement-error model. They do not establish that a
grid is an appropriate scientific scoring resolution.

Before this supplement, genuine arithmetic defects were corrected and prior
attempts retained:

1. Exact constant decimal vectors are recognized before variance/correlation.
2. Averaging identical predictions over one/two/three alleles preserves the
   constant exactly.
3. Source observed-position ranks use exact rational decimal sums/counts.
4. Identical raw feature rows reuse an identical computed prediction within a
   fit.

These reduce spurious differences but do not force distinct-input predicted
near-ties to match. No epsilon was inserted into frozen v1 exact ranking.

## Artifact audit

The script
[analyze_rank_sensitivity.py](../scripts/analyze_rank_sensitivity.py):

- Verifies the frozen supplement and parent protocol hashes.
- Verifies independent artifact and original cohort commitments.
- Matches all 131112 prediction identities and their source observations.
- Recomputes all 480 native archived correlations before evaluating grids.
- Retains all 1440 metric-pair/level rows, every undefined status, all differences
  and per-level worst cases.
- Performs no network access, model fit or selection.

The original diagnostic ran at `/tmp/helix-rank-sensitivity-v1`. Its artifact-index SHA-256 is
`b8a43b694eba11f1a95a93c0b6574e0415d9d31894a7770201d94caef279daa1`.
That original bundle is archived in the published snapshot's `history/`;
the clean-revision re-execution is in
[`results/regulatory-v1/rank-sensitivity/`](results/regulatory-v1/rank-sensitivity/summary.json).
It includes `summary.json`, `all-metric-pairs.json`,
`raw-archive-checks.json` and hash commitments. Complete provenance records:

| Input | SHA-256 |
| --- | --- |
| Rust exact-source study | `caf376dcb069893767f0b010afc2b5b56a7099f657890316a2977dfdbe34e507` |
| Independent artifact index | `83e8f938260bbcb4f3dd73912eaee1637a3855aff9021a3f9eac17abea20ea88` |
| Cohort manifest | `200342ba098e4e302cb306af67f68b01fccf13291ad6650b7dcb1281d83ddae5` |
| Diagnostic implementation | `31bc0e81ac0bfb9d440ce28021c77dfd1f4d768952efdb21ff454951413e013e` |

The diagnostic ran with Python 3.12.14 and NumPy 2.3.5. The stored script accepts
plain or gzip-compressed Rust JSON and refuses an existing output directory.
Synthetic controls check exact positive/negative half-grid ties, source-decimal
ties and undefined constant-vector ranks.

## Consequence for the research report

Lead with the independently reproduced primary loss results and their negative
predictive conclusion. Preserve raw secondary-rank failures. The finite-grid
tables may explain numerical sensitivity, but they cannot rescue a model that
did not improve the primary error, establish a mechanism, or count as genetic
novelty. This is an auditable methods finding and a reason to predeclare rank
reporting precision in a future study.
