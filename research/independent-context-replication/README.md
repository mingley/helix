# Independent context-specificity replay

The frozen `context-specificity-v2` replay passes its full independent numerical
comparison. The source values do not obey the cited simple zero-sum ESS formula;
this source-version/normalization issue remains unresolved. Values were retained
unchanged. Read `review-v2.md` for the result and limits.

The implementation uses CPython 3.12.14 and NumPy 2.3.5 from the pinned research
environment. It imports neither Rust code nor the previous regulatory replay.
It independently parses raw XLSX, validates source eligibility and every
declared-strand hg38 interval, reconstructs families and guards, extracts
features, solves ridge fits and replays resampling.

## Commands

Run from the repository root after installing the pinned research environment
and obtaining the checksum-verified reference described in
`research/context-reference/README.md`. Supply a group ledger bound to this panel
and protocol. Every output directory must be new; previous attempts are retained.

```sh
source scripts/research-env.sh
OPENBLAS_NUM_THREADS=1 python3 -m unittest discover \
  -s research/independent-context-replication \
  -p test_context_replication.py -v

OPENBLAS_NUM_THREADS=1 python3 scripts/replicate_context_study.py \
  --protocol experiments/context-specificity-v2.json \
  --dataset datasets/context-panel \
  --groups runs/context-integration-v2-groups.json \
  --reference /workspace/.tools/reference/hg38.2bit \
  --out runs/context-replay-example-python

OPENBLAS_NUM_THREADS=1 python3 scripts/replicate_context_study.py \
  --compare-rust runs/context-integration-v2-rust-corrected/study.json.gz \
  --independent-run runs/context-replay-example-python \
  --out runs/context-replay-example-comparison

python3 scripts/replicate_context_study.py \
  --audit-endpoint \
  --independent-run runs/context-replay-example-python \
  --out runs/context-replay-example-endpoint-audit
```

The comparator exits 0 on a complete match, 2 on numerical/membership mismatch,
and nonzero on malformed input or provenance error. It detects gzip by magic
bytes and records compressed and decoded commitments. The current schema-4 core
must contain the source-feature dump; the initial core predating it remains a
historical attempt.

`--source-check-only --out NEW_DIRECTORY` verifies source data without fitting.
`--describe` prints endpoint, protocol hash and runtime versions. Outputs are
never overwritten.

## What is independently checked

The only accepted protocol SHA256 is
`8f0576845fecb5cfcf9bfa354fab446d849244448a3489905a6ce33bdeb21cec`.
Declared file hashes and semantic gates are checked before parsing outcomes.
ZIP/XML cell types and literals are read without the root importer or openpyxl.
Strings, rich text, numeric literals, booleans and errors remain distinct;
formulas fail.

All 60,000 S7 design IDs and 55,338 S9 score IDs are considered. Exact joins,
finite complete triples, source230[15:215] core extraction and declared-strand
hg38 eligibility reproduce all 55,276 admitted elements. Missingness and
exclusions are retained without imputation.

Families use bytewise NumPy Hamming verification instead of root's packed
integers. All 372,850 source-ID candidates, unions and 2,761 component hashes
match. Guards search merged expanded test intervals instead of root's
prefix-maximum algorithm. Every five-fold train/test/whole-component guard list
matches.

Features are integer RC-canonical 3-mer counts divided once by 198, then GC count
divided once by 200. All integer counts and 1,824,108 feature fractions match the
Rust dump exactly. Features use the core only.

Ridge uses training-only population SD, unit scale for SD at most 1e-12, free
intercepts and alpha 10 summed squared error. Shared slopes fit the three-context
mean response per training element; independent context fits use the same
sequence design and alpha. Additive and shared models reuse the same slopes.
No response rescaling, rank endpoint, tuning, feature selection or zero-sum
correction occurs.

All five models, three fixed comparisons, context pairs and 2,000 common
SplitMix64 draws are retained. Primary resampling pools loss sums and element
counts with component multiplicity rather than averaging component means.
Equal-component weighting is a fixed sensitivity.

The comparator binds core inventory and source/protocol/panel/group commitments.
It checks observations, memberships, fitted statistics, features, predictions,
all sequence/group/context/pair losses, the common draw table, all-model sampled
losses and percentile intervals. Prediction/fit tolerance is 1e-7 and
loss/statistic tolerance is 1e-8; observations and feature fractions match
exactly. No tolerance was widened.

Artifact hash verification does not replace the CLI ledger audit. The initial
Rust core had a separate serialization-order integrity-audit failure; the
corrected core retains identical targets/objectives and passes this comparison.

## Tests and retained attempts

All 34 synthetic tests pass: hand-counted features, closed-form ridge,
unpenalized-context normal equations, training-only scaling, opposite slopes,
aligned penalties, balanced identities, test-target perturbation, missingness,
exhaustive small-graph grouping, inclusive 10kb guards, candidate-budget
accounting, XLSX types, hexadecimal seed, unequal-size component sampling, exact
source-sum arithmetic, gzip commitments and artifact mutation/path failures.

`run-summary-v2.json` records input/output hashes and all retained attempts: one
successful Python fit, initial/corrected Rust cores, both comparisons and the
separate exact-rational endpoint audit. Full compressed outputs remain in their
run directories for packaging; compact comparison and endpoint-audit reports
are retained here.
