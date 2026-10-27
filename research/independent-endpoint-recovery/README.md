# Independent endpoint-recovery verifier

The independent openpyxl implementation reproduced the primary source audit:
**7,008 checks passed** for the original 36 routes and direct ESS mapping;
**1,162 checks passed** for the separate, explicitly post-hoc six-route audit.
All 36 original activity-to-ESS routes and all six later routes failed their
frozen source-agreement tests. Published S11 ESS reproduced S12 after sample
standardization on all 55,338 exact S12 IDs, with maximum independent error
6.217248937900877e-15 against the frozen 1e-12 limit.

This establishes numerical compatibility without identifying producer history.
The tool performs no predictive fit, acquisition, imputation, target replacement
or selection of a best route. Previous results remain unchanged. The complete
[published evidence](../results/context-endpoint-recovery-v1/README.md) retains
the arrays, source ledgers, implementation snapshots and failed comparison
attempts. See [review.md](review.md) for commitments and scientific limits.

The original protocol SHA256 is
`00175e41cf9fa0cbf177ce8fc7ef143eb508aec8b11d058f7acb4458284dc868`;
the separate six-route protocol SHA256 is
`a0e3f521d1086fca0633fc4c473a3655aaf2fa9d4a4818b0429a735f46d42372`.
Both were frozen before their respective computations. The latter protocol was
written after the original results and cohort gap were known, and remains a
post-hoc diagnostic.

Exact source/header binding excludes unbound numeric fields. Worksheet
dimensions are advisory; physical stored rows, raw numeric literals, formulas,
types, missingness and exact case-sensitive IDs are preserved. All eight bound
sheets and 796,645 stored data rows matched the primary strict-XML export.
Formula caches, booleans, errors and blanks cannot silently become measurements.
Every route retains its declared normalization population; no intersection,
imputation or missing-value mean rescues a failed or non-estimable route.

Use the pinned research environment and fresh output directories for every
replay. The example paths below must not exist beforehand. Run from the
repository root; the published primary audit is read without alteration.

```sh
source scripts/research-env.sh
export PYTHONDONTWRITEBYTECODE=1 OPENBLAS_NUM_THREADS=1
python3 scripts/verify_context_endpoint_recovery.py \
  --empirical \
  --primary-export research/context-endpoint-recovery/source/jointdata-source-values.json.gz \
  --out runs/endpoint-independent-replay
python3 scripts/verify_context_endpoint_recovery.py \
  --compare-primary research/results/context-endpoint-recovery-v1/primary \
  --independent-run runs/endpoint-independent-replay \
  --out runs/endpoint-comparison-replay
python3 scripts/verify_context_endpoint_recovery.py \
  --ess-cohort-audit \
  --independent-run runs/endpoint-independent-replay \
  --recovery-comparison runs/endpoint-comparison-replay/comparison.json \
  --out runs/ess-cohort-independent-replay
python3 scripts/verify_context_endpoint_recovery.py \
  --compare-posthoc research/results/context-endpoint-recovery-v1/addenda/context-ess-cohort-audit-v1 \
  --independent-run runs/ess-cohort-independent-replay \
  --out runs/ess-cohort-comparison-replay
python3 -m unittest discover \
  -s research/independent-endpoint-recovery -p test_endpoint_recovery.py -v
```

All **28 synthetic controls passed**. They cover physical-row parsing, raw
literals and typed cells, source commitments, exact joins, tolerance boundaries,
population and replicate arithmetic, all route tuples, undefined variants,
non-identifiability and exact ESS-cohort guards.

The pinned reader uses CPython 3.12.14, NumPy 2.3.5, openpyxl 3.1.5 and
et-xmlfile 2.0.0. Official PyPI wheel hashes, byte sizes and URLs are retained in
[dependency-receipts.json](dependency-receipts.json); the central research lock
pins installation. Metadata-only inventory is available with `--xlsx PATH
--expected-sha256 HEX --header 'Sheet:1' --out NEW`. An optional MD5 must be
independently declared; an opaque publisher ETag is never interpreted as MD5.
Hashes establish consistency with commitments, not publisher-signed authenticity.
