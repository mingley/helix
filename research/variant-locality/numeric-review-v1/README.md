# Independent numerical review of frozen variant predictions

An independent standard-library implementation parsed the published activity-v3
fits and reproduced **24,960 ensemble predictions** without importing the
predictor, fitting a model or reading new experimental outcomes. The comparison
used the exact rational interpretation of every saved binary64 coefficient and
training scale. Maximum ensemble difference was **5.551115123125783e-17**;
maximum difference from direct subtraction of two standardized full-sequence
predictions was **8.066464163292153e-16**.

There are 4,992 synthetic SNVs per predictor configuration: all 3,072 interior
5mer/alternative-base combinations and all 1,920 boundary cases. The latter
enumerate 192 cases each at positions 0 and 199 and 768 each at positions 1 and
198, using zero-based coordinates. The five configurations are the independent
cell fits `I/HepG2`, `I/K562`, `I/WTC11`, the shared `S` fit, and additive `A`.
`S` and `A` use the same DNA slopes and give identical substitution deltas;
their context labels in this review are both K562. These are five predictor
configurations, not five different feature architectures.

All **24,960 exact ensemble comparisons** pass the review's absolute `1e-12`
limit. The review also records the maximum discrepancy across all **124,800
fold-level comparisons**, which is `5.551115123125783e-17`. The predictor's
separate direct-subtraction guard is `1e-10`; all predictions pass it. Exact
zero feature changes produce zero predictions in both implementations.

| Predictor | Maximum absolute interior prediction | Maximum absolute boundary prediction |
| --- | ---: | ---: |
| I/HepG2 | 0.10294832690298229 | 0.125604494439756 |
| I/K562 | 0.13234503481228885 | 0.15373186010855194 |
| I/WTC11 | 0.1288579581406721 | 0.150215233589851 |
| S | 0.09159806904071756 | 0.13250841807448463 |
| A | 0.09159806904071756 | 0.13250841807448463 |

These are exhaustive bounds on these particular frozen ensemble predictions
for same-length 200nt SNVs, expressed on the transferred published
library-normalized log2(RNA/DNA) model scale. The ranges are symmetric under
allele reversal. Boundary predictions can exceed the interior bounds, so an
interior-only envelope must not be described as a bound for every position.
The full exact rational extrema are retained in [summary.json](summary.json).
No observed effect sizes, calibration, biological limits or clinical conclusions
follow from these model bounds.

Each configuration has 1,405 distinct interior prediction values, including 64
zero edit cases. Each has 896 distinct boundary prediction values and no zero
boundary cases in this enumeration. These are **counts of predicted values**
under the specified coefficients; this review does not assert that 896 is a
count of distinct boundary feature vectors.

For a fixed fold and sequence length, the exact mathematical bridge is
`prediction(x) = B + sum_j W_j count_j(x)`, where
`W_j = beta_j / (training_scale_j * denominator_j)`,
`B = intercept - sum_j beta_j * training_mean_j / training_scale_j`, and the
denominators are 198 for 3mer frequencies and 200 for GC fraction. A five-fold
ensemble averages these weights. Finite binary64 inputs are exact rationals,
so a common positive denominator can convert that mathematical score to an
integer-weight score. The Lean locality proof applies to that exact
architecture. It does not prove identities about IEEE floating-point execution:
division, products, addition and full-prediction subtraction still round.
This review measures their discrepancy separately. The predictor computes
count differences before division, so identical count-delta vectors follow
identical arithmetic regardless of unchanged remote sequence context.

The original `/tmp` computation is retained byte-for-byte: `implementation.py`,
`summary.json`, all three pair-input files, all three primary output folders and
their original `manifest.json`. That original manifest binds the original
eleven evidence files; its SHA-256 is
`b700ae1057746416e6559d65dd97dc4551e523463fefcde4cc0129fe2e3f23ab`.
`reviewed-predictor.py` is an additional byte-exact snapshot of the reviewed
predictor. `publication.json` inventories every file in this published folder
except itself, including this documentation and both implementations.

The frozen training artifact manifest is
`2ea476bc5413dc13486dec9f3d7e9162a50a0e47f1d83ce29efabc5f5cc5af9d`.
The reviewed predictor SHA-256 is
`7e6123eab6c8c55cd714bd18b092bb8e5806f87a1640a44869404ae8dd30fae4`.
The independent script verifies every artifact commitment before parsing the
fits. All inputs here are synthetic DNA and already-published training fits;
no new external experimental response payload was acquired or inspected.

Replay from the repository root with a fresh output directory and the reviewed
predictor still at `scripts/predict_variant_effects.py`:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 research/variant-locality/numeric-review-v1/implementation.py \
  --root "$PWD" --out runs/variant-numeric-independent-replay
```

The script imports no predictor code; it runs the predictor CLI as the object
under review, then independently computes exact rational predictions. A future
change to that CLI is a new review target and may change outputs. The retained
`reviewed-predictor.py` can also be run directly with each retained `*-pairs.json`
file, model `I`, `S` or `A`, the frozen core and manifest pin, and a fresh
`--out` directory. All replay output paths must be new.
