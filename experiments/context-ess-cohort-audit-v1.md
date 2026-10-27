# Post-hoc audit of the published ESS cohort

The frozen authority is
[context-ess-cohort-audit-v1.json](context-ess-cohort-audit-v1.json), SHA-256
`a0e3f521d1086fca0633fc4c473a3655aaf2fa9d4a4818b0429a735f46d42372`.
The parent committed it as `bc766a2` before any activity normalization on the
exact S11 ESS cohort. This is a separate post-hoc arithmetic diagnostic, not
a predictive study or a revision of the original 36-route results.

## Exposure and question

All earlier source values, models, arithmetic findings and 36 reconstruction
results are known. The published activity aggregate has 56977 complete
three-cell rows, whereas the final S11 ESS table has 55977. The first audit
normalized the activity bases on their own prescribed source populations;
none of its 36 variants reconstructed published S11 ESS within tolerance.

This audit asks one narrower question: does normalization on **every exact S11
ESS ID**, rather than the larger activity population, produce a compatible
reconstruction? The target values are already known. There is no preregistered
or independent-validation claim, and no requirement to obtain a favorable match.

## Fixed inputs and six routes

Use the three unchanged, hashed activity arrays from the completed primary
recovery audit:

1. Published activity aggregate.
2. Arithmetic mean of three published replicate scores.
3. Arithmetic mean after explicit per-replicate median subtraction, preserving
   the medians originally computed over all finite source values per replicate.

Do not recompute replicate medians on the new cohort or change empirical values.
Normalize each base on the complete, unique, finite S11 ESS population of
**55977 IDs**, ID-set SHA-256
`3f5667fe22887a579d66da4b0b4e2e59f2196f0d6d465d0656a10ef2b6fe3c8f`.

Each activity base must cover **every** target ID with three finite values.
If it does not, both routes for that base are nonestimable. Report missing IDs
and contexts; do not use an intersection, partial mean or imputation.

For each eligible base, compute per-cell z-scores on this fixed population using
`ddof=0` and `ddof=1`, then subtract each element's three-cell mean. These are
the only six routes. An exactly constant context, nonpositive/nonfinite SD or
insufficient population makes a route nonestimable. Exact all-equal checks
prevent tiny floating SDs from masquerading as variability.

Compare every reconstructed array to unchanged S11 ESS on the full cohort,
the exact old-S9 intersection and the fixed old admitted intersection. The
general rule remains

`abs(candidate-source) <= 1e-8 + 1e-8*abs(source)`.

Report all parameters, counts, ID hashes, signed errors, MAE, RMSE, maximum
error, Type-7 quantiles, within-tolerance counts, column means and row sums.
Derived floating arrays have no source decimal lexeme: exact Decimal equality
is not claimed for them. Preserve all six results and all nonestimable reasons.

Because every reconstructed column was centered on this same population, its
mean should be approximately zero. Report the corresponding published S11
means without forcing them to zero. This structural check can reject a proposed
identity independently of whether a fitted-looking example appears close.

## Fixed S12 derivative

For **every estimable route**, whether or not it matches S11, restrict to every
exact S12 ID: 55338 source IDs, with full coverage required. Subtract each
column mean and divide by sample SD (`ddof=1`) on that full S12 subset.
Compare to unchanged S12/old S9 using maximum absolute tolerance `1e-12`, on
all S12 IDs and the fixed admitted intersection. Missing coverage blocks the
derivative; no favorable subset is allowed.

This sample-z-score step is the one previously specified release-comparison
operation. It is not a freely fitted affine calibration. Do not add weights,
offset searches or additional transforms. Some routes can become equivalent
after this step; report equivalence rather than selecting a unique history.

## Metadata and interpretation

Publish full ID differences/intersections and counts by the existing exact S10
`category` field, including unmapped IDs. Categories describe the cohort gap;
they cannot choose, remove or weight the normalization cohort. Do not infer
new activity-based categories.

A compatible arithmetic route does not prove the producer used it. Per-cell
offsets and positive scale changes can disappear under z-scoring; replicate
median offsets and some SD conventions are observationally indistinguishable.
Keep historical identification, source-file identity and numerical agreement
separate. Earlier 36-route negatives, old targets, fits and losses remain
unchanged. No predictive fitting, independent biological validation or regulatory
novelty claim is authorized.

The implementation is
[audit_context_ess_cohort.py](../scripts/audit_context_ess_cohort.py). It verifies
all frozen inputs and refuses existing output directories. For portable replay,
`--primary-root` explicitly supplies a different directory containing the same
frozen primary files by basename; every byte count and SHA-256 still must match.
`--source-workbook` similarly supplies the byte-identical official workbook.

## Completed diagnostic

All six routes were estimable with every required target ID and all six failed
both reconstruction comparisons. No empirical values, earlier routes or
predictive fits were changed.

| Activity base | Pre-SD ddof | Maximum error versus S11 | Maximum error versus S12 derivative |
| --- | ---: | ---: | ---: |
| Published aggregate | 0 | 1.57640265 | 1.20917329 |
| Published aggregate | 1 | 1.57638057 | 1.20917329 |
| Three-replicate mean | 0 | 1.57690317 | 1.20917261 |
| Three-replicate mean | 1 | 1.57688108 | 1.20917261 |
| Median-subtracted replicate mean | 0 | 1.57690317 | 1.20917261 |
| Median-subtracted replicate mean | 1 | 1.57688108 | 1.20917261 |

None of the 167931 S11 comparison cells satisfied the frozen general tolerance
for any route. All six S12 derivatives failed the maximum-error `1e-12` test.
The separately tested transformation of **published S11 ESS itself** to S12
remains verified by the earlier recovery audit; the unresolved link is from
the tested activity bases to S11 ESS.

The 1000 IDs present in the complete published activity cohort but absent from
S11 ESS span existing source categories: 193 HepG2, 323 K562 and 322 WTC11
potential enhancers; 118 promoters; 33 synthetic positive controls; six shuffled
negative controls; four synthetic negative controls; and one K562 positive
control. No category was used to change a cohort. Neither this cohort difference
nor these six formula variants recover the missing upstream processing history.

The full [summary](../research/results/context-endpoint-recovery-v1/addenda/context-ess-cohort-audit-v1/summary.json) SHA-256 is
`74083f00da45de823bb2f10177bbe7c7c010e0f14a5291db32d03396b9d18906`;
implementation SHA-256 is
`f2267f24a6f3c0f3d2b6bb7ec23865df74d1e63ede0b12766fd4d435e586d335`.
All 20 artifact commitments verified. Five synthetic controls passed, and an
explicit transport-root replay produced a byte-identical summary. Independent
numerical comparison passes all 1,162 checks, reproducing all six negative
results. Its [receipt](../research/results/context-endpoint-recovery-v1/addenda/context-ess-cohort-audit-v1-comparison-v2/comparison.json)
records the explicit reserved missing-category label adapter; the earlier
metadata-only failed comparison is retained in publication history.
