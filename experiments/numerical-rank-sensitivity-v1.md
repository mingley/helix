# Post-hoc numerical rank sensitivity v1

This is a separately versioned methods diagnostic. It does not change the frozen
regulatory prediction experiment or turn its failed exact-rank endpoint into a
pass. The executable supplement is
[numerical-rank-sensitivity-v1.json](numerical-rank-sensitivity-v1.json), SHA-256
`b968f0b3519b3890c36d2756ce97780f3a42effb8ee06b138e4cf395138b89ed`.
The complete grid was fixed before this diagnostic was run.

## Why the supplement exists

The two implementations agree on predictions to much better than the frozen
`1e-7` absolute tolerance, and on primary errors, splits, fits and controls.
However, Spearman is discontinuous at ties. Tiny floating-point differences can
split mathematically tied predictions and change ranks even when they have no
meaningful effect on squared error. Source-decimal reductions can also create or
break artificial position-mean ties if calculated naively.

The native v1 exact-rank metrics remain recorded with their exact failures.
Exact source-decimal arithmetic and identical-feature prediction caching are
separate semantic corrections already documented. This supplement asks a new
question: how do the reported rank comparisons depend on explicitly declared
prediction-reporting resolution?

## Fixed estimands

Report every one of these levels:

| Level | Prediction reporting quantum, log2-effect units |
| --- | ---: |
| raw exact | No quantization; original saved floating values |
| 1 | `1e-12` |
| 2 | `1e-10` |
| 3 | `1e-8` |
| 4 | `1e-6` |
| 5 | `0.01` |

The grid spans fine arithmetic differences through the source's reported
coefficient scale and brackets the existing prediction-reproduction tolerance.
No level is selected according to agreement or apparent biological association.

Quantize **predictions only**, after each implementation's native row/position
reduction. Convert the saved binary float to an exact rational, divide by the
exact decimal quantum, round to the nearest integer with ties to even, and rank
the integer codes directly. This avoids adding another binary-division rounding
ambiguity. Do not convert codes back to floating grid values before ranking.

Observed row ranks use original source decimal scores. Observed position means
are exact rational sums of those decimals divided by available-allele count;
they are never quantized. At the position level, Rust uses its saved native
position predictions, while Python reproduces the archived constant-preserving,
ALT-sorted NumPy reduction. The script first checks that these reconstructed raw
metrics agree with each implementation's own archive.

Average exact ties and return an undefined correlation if a rank vector is
constant. Quantization can increase ties or make a metric undefined. Such loss of
resolution is reported explicitly, not hidden inside an agreement percentage.

## Complete reporting

Compare all 240 paired metrics: four evaluation families, six models and every
held-out assay, with row and position correlations separately. For each of the
six levels retain:

- Both correlations, defined/undefined status and absolute difference.
- The number of targets and distinct observed/predicted rank inputs.
- The maximum original prediction difference and both original prediction
  ranges.
- Every discrepancy, not just the maximum.
- Counts, maximum differences and ten worst identities at each level.

The inherited `1e-8` metric tolerance is only a reporting comparator. No new
p-value, biological significance threshold, model fit or feature selection is
introduced. Coarse-resolution agreement is a different estimand, not retroactive
success for raw exact Spearman.

Input artifact hashes, original source hashes, parent/supplement hashes and
software/code versions are retained. Existing outputs are never overwritten.

## Run

```bash
source scripts/research-env.sh
python scripts/analyze_rank_sensitivity.py --self-test
python scripts/analyze_rank_sensitivity.py \
  --rust /path/to/final-rust-study.json.gz \
  --independent /path/to/final-independent-directory \
  --out runs/rank-sensitivity-new
```

An uncompressed Rust JSON is also accepted. This script performs no network
access and no model fitting. The result is a numerical-methods finding, not
genetic novelty or a mechanism explaining reporter biology.
