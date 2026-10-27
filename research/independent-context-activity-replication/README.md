# Independent activity-v3 replay

The frozen activity replay passes **258/258 numerical comparison categories**,
including exact published source observations, **1,824,108 exact feature
fractions**, all fitted quantities, **829,140 model predictions** and every
declared loss/resampling/common-activity check. Maximum prediction difference
is **6.375e-12**, below the frozen absolute 1e-7 tolerance. Metrics use absolute
1e-8 tolerance; neither was widened. See the
[activity study](../context-activity-study.md) for the scientific result.

The primary C-minus-A gain is 0.064171577, approximately 8.7996% of the
context-mean MSE, on the mean of three published library-normalized activity
scales. It is not absolute expression, universal TF strength or independent
biological replication. Upstream activity-to-ESS history remains unresolved.

## Source gates and independent computation

Authority: `experiments/context-activity-v3.json`, SHA-256
`64e14037fb1f42ddc1ecaca1d23bb8096fabac1f1e892d105eaa20f525fbbb20`.
The protocol was frozen after archive acquisition began but before inspecting
new workbook headers/activity values. Prior specificity results, public methods,
code/schema expectations and panel selection were known.

Before opening the activity workbook, the adapter verifies:

- The dataset manifest, protocol/endpoint, every committed receipt and all five
  source gates: archive integrity, empirical schema, complete coverage, physical
  design identity and recorded aggregation uncertainty.
- The original S11 workbook's 53,557,805 bytes and SHA-256
  `9c95ca2b8df0ff3edc7dda4fdcf337f33e652b96a0833a70e786bbcdc2b76b4e`.
- Exact prior DNA-panel/grouping commitments and all 55,276 fixed IDs. Old ESS
  response fields are discarded; DNA features are extracted anew.

Pinned openpyxl 3.1.5 independently reads the raw `all_cell_types_summary`
worksheet, binding `name` and the literal three `<cell> [log2(rna/dna)]`
columns. Formula, error, boolean, numeric-text, missing or nonfinite required
values block the replay. IDs are exact; there is no suffix/fuzzy rescue,
subset selection, imputation, ESS inversion or outcome rescaling. The root
normalized activity export is hashed as provenance but is not the target source.

The replay reuses inspected independent NumPy fitting/bootstrap functions from
the earlier Python implementation, not Rust calculations. All new outcome fits,
training-only scalers, errors and resamples are recomputed. Reused metadata
comprises the previously independently validated 200-base cores and frozen
2,761 groups/five guarded folds; v3 does not claim another independent reference
or group reconstruction. Official final design-byte identity supports that reuse.

All five models and three comparisons remain in the original array order:
S-minus-A diagnostic, C-minus-A primary at index one, A-minus-I secondary.
Primary loss weights elements equally after averaging their three cells;
bootstrap draws resample components and pool member loss sums/counts. A generic
Rust-core note calling group weighting “primary” is inapplicable to v3; the
frozen caller protocol and reported computations use **equal-element weighting**.
Equal-component metrics remain a sensitivity.

## Commands

After the pinned setup, run synthetic controls without source outcomes:

```bash
source scripts/research-env.sh
OPENBLAS_NUM_THREADS=1 python3 -m unittest discover \
  -s research/independent-context-activity-replication \
  -p test_context_activity_replication.py -v
python3 scripts/replicate_context_activity.py --describe
```

The full workflow verifies source gates, runs Rust/audit, independent Python,
comparison and three figures. Choose a new prefix:

```bash
bash scripts/check_context_activity.sh runs/activity-review
```

Or run the independent replay/comparison explicitly against an existing Rust
core, with new output directories:

```bash
OPENBLAS_NUM_THREADS=1 python3 scripts/replicate_context_activity.py \
  --out runs/activity-review-python
OPENBLAS_NUM_THREADS=1 python3 scripts/replicate_context_activity.py \
  --compare-rust runs/activity-review-rust/study.json.gz \
  --independent-run runs/activity-review-python \
  --out runs/activity-review-comparison
```

`--source-check-only --out NEW_DIRECTORY` verifies/parses source without fitting.
It still requires the passed dataset gate. The comparator exits 0 for a full
match, 2 for a comparison mismatch, and nonzero for malformed/provenance input.
Output directories are never overwritten. Offline publication inspection is:

```bash
bash scripts/check_context_activity.sh --verify-published
```

## Controls, artifacts and retained failure

All **17 new synthetic controls** pass. They test exact headers/IDs and numeric
types, every pre-source gate, mandatory uncertainty recording, additional
committed provenance, missingness/duplicate rejection, a hand-computed ridge
gain, loss of the common component in a centered toy example, and primary
bootstrap ordering. No empirical workbook is used by these tests.

The initial `runs/context-activity-v3-final-python` attempt stopped **before raw
activity parsing or fitting** because the verifier incorrectly required exactly
seven receipt filenames. The finalized manifest had additional valid provenance
files. The corrected verifier requires essential receipts, accepts safe additional
names and checks every listed hash. The failed attempt is retained; targets,
models, penalty, folds, metric definitions and tolerances were unchanged.

Successful replay: `runs/context-activity-v3-final-python-corrected`.
Comparison: `runs/context-activity-v3-final-comparison`.
The [published snapshot](../results/context-activity-v3/README.md) retains the
successful directory under `independent-python/`, with failed history separately
recorded.

Each manifest binds source/protocol/group commitments, parser and model-code
snapshots, runtime/dependency hashes, exact source records, all fit/prediction
arrays, common-activity diagnostics and the common bootstrap table. Integrity
hashes establish byte consistency, not source authenticity or biological truth.
The common-activity and pairwise identities are correlated implementation
checks, not extra biological experiments. All conditional ranges omit refitting,
source-normalization uncertainty and study-to-study variation.
