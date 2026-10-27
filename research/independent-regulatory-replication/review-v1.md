# Independent replication review: frozen v1

## Result and scope

The independent NumPy/SciPy implementation verified the primary model
calculations and numerical controls. It did **not** reproduce every secondary
exact-rank statistic. The full comparison remains failed, with an exit status
of 2; this is not relabeled a complete pass.

| Check | Outcome |
| --- | --- |
| Feature vocabularies | All six models agree |
| Predictions | 131,112 target predictions agree within frozen 1e-7 absolute tolerance |
| Maximum prediction difference | 7.729e-12 |
| Train/test/guard memberships | All 84 lists agree |
| Fit preprocessing/coefficients | All 168 fits agree within frozen metric tolerance |
| Macro MSE/MAE | All 24 model/family results agree, maximum MSE difference 3.03e-15 |
| Per-assay losses | All 120 model/assay results agree |
| Bootstrap draws | Both 2,000-draw tables match exactly |
| Bootstrap statistics | Maximum difference 9.04e-15 |
| Spatial-null contrasts | All six sets of 99 contrasts agree, maximum difference 1.39e-14 |
| Exact secondary Spearman | 26 failures remain: 24 row-level and two position-level |

See `comparison-v1.json` for every discrepancy and `run-summary-v1.json` for
input/software commitments, frozen comparisons, context summaries and the
retained attempt inventory.

This establishes agreement of independent software on shared inputs, not
independent biological replication. The nine primary elements share a study,
assay design and data processing. HBB had already been explored.

## What the primary results support

Macro MSE is position-weighted within each assay and equal-weighted across
primary assays. Units are squared log2 reporter-effect units.

| Model | HBB guarded CV | Nine-promoter leave-element-out | HBB row-random diagnostic |
| --- | ---: | ---: | ---: |
| Mean | 0.021933886 | 0.116392187 | 0.021462133 |
| Chemistry | 0.022683072 | 0.116899210 | 0.021725836 |
| Coordinate | 0.026119801 | 0.116528763 | 0.020792210 |
| Delta 3-mer | 0.025027226 | 0.118370897 | 0.020941712 |
| Motif | 0.024996569 | 0.116829487 | 0.021855343 |
| Full | 0.029609654 | 0.117988875 | 0.020436205 |

No sequence model beats the mean baseline on the primary cross-element macro
error. The full model's lower row-random HBB error does not survive the guarded
within-element test. The split difference also changes difficulty and cannot
be attributed solely to leakage.

The motif-versus-chemistry primary comparison has a small positive improvement,
6.97e-5, with a descriptive element-bootstrap interval
[-5.80e-4, 6.43e-4]. Its chemistry-stratified spatial-null tail is 0.01,
Holm-within-family 0.03, and conservative six-tail Bonferroni 0.06. These are
frozen robustness diagnostics, not calibrated biological-significance values.
The motif model still loses to the mean baseline; this does not establish a
new regulatory mechanism or superiority to standard genomic predictors.

## Why the exact-rank check remains failed

Initial attempts exposed real numerical defects, which were corrected without
changing models or tuning to outcomes:

- Constant predictors could acquire artificial position ranks when one-, two-
  and three-allele groups were re-averaged.
- BLAS batch alignment could give identical feature rows predictions differing
  by one ULP.
- Floating regrouping could split exact source-decimal means before ranking.

The final implementation preserves exact constants, evaluates each identical
unscaled feature row once, and ranks observed position means using exact
rational arithmetic from raw decimal scores.

Remaining failures are consistent with exact-rank discontinuity near
mathematically tied predictions from **different** feature vectors, plus the
declared floating reduction conventions. A concrete case is the HBB coordinate
fit: A>G and A>T each contain 36 observations with the identical exact score sum
2/5 and mean 1/90. The ridge problem has a corresponding feature-exchange
symmetry. Independent floating solvers can split the equal coefficients by
approximately 5e-18. Their predictions agree to far better than the frozen
tolerance, while exact rank order can differ.

No coefficient rounding, rank epsilon, score grid or relaxed tolerance was used
to force agreement. A separate post-hoc sensitivity analysis may describe this
instability; it cannot retroactively change the v1 pass criterion.

## Orientation and biological limits

Motif calculations scan both DNA strands and pass reverse-complement invariance
checks. HBB's physical assay insert is native-minus relative to the source's
plus-strand coordinate target. That orientation distinction does not change
these motif deltas. It does not establish the reporter's actual transcription
start or RNA sequence.

The fragment includes native 5-prime UTR sequence and upstream sequence, so
neither the source name “promoter” nor negative coding coordinates assign all
bases a promoter-only mechanism. Reporter effects are not endogenous HBB
expression, hemoglobin function, or clinical pathogenicity.

## Reproduction and retained attempts

The final independent directory is:

    /tmp/helix-independent-regulatory-v1-exact-source

Its exact Rust counterpart is:

    /tmp/helix-regulatory-rust-v1-exact-source.json

The final comparison is:

    /tmp/helix-regulatory-independent-comparison-v1-exact-source/comparison.json

The immutable output inventory and reasons for earlier attempts are recorded in
`run-summary-v1.json`. Root archives complete result bundles separately.
Temporary runtime paths identify the actual executed attempts; new users should
run the documented commands with new output paths.

The final Python implementation has 18 synthetic golden/regression tests.
They verify model arithmetic, sparse-position eligibility, feature counts,
guard geometry, PRNG streams, rotation order, dual-strand motif scoring and
the numerical corrections. None requires fitting real assay outcomes.
