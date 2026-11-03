# Post-hoc feature collisions in the published activity panel

This exploratory audit reuses **only the already-analyzed v3 reporter-activity
labels** and the existing DNA panel. It finds **no equal-feature pair of
distinct, non-reverse-complement 200-base cores** among 55,276 panel records.
There is exactly one feature-collision class: two distinct cores that are whole
reverse complements. All 55,276 complete 33-integer feature vectors, and their
normalized floating values, were rebuilt from source DNA and match the retained
v3 feature snapshot exactly.

The two colliding records have different published aggregate activity labels in
all three measured cells. The exact equal-feature lower bound on their
**observed-label** mean squared error is approximately **0.040247** averaged
over those two records and the three cells. When every one of the 55,276 panel
elements and three cells has equal weight, the corresponding entire-panel floor
is approximately **0.00000145622**. These are finite-sample arithmetic bounds,
not estimates of irreducible biological variation, new laboratory replication,
or held-out generalization. Measurement noise, prior source normalization and
this single pair's shared source context remain part of the observed labels.

## Cohort and exact equivalence

The source is the previously published
[context-activity-v3 study](../context-activity-study.md), which analyzed
55,276 source-selected 200-base cores across HepG2, K562 and WTC11. Its target
is the unchanged published library-normalized log2(RNA/DNA) aggregate from
Agarwal et al., *Nature* (2025),
[doi:10.1038/s41586-024-08430-9](https://doi.org/10.1038/s41586-024-08430-9).
Those 165,828 labels have already been analyzed in that study; this audit was
designed **after** that exposure. It does not read Siraj or Kreevan outcomes or
claim a new outcome-naive experiment. The source's upstream normalization and
precise historical replicate aggregation remain subject to the qualifications
in the v3 report.

Features are the 32 reverse-complement-canonical overlapping 3-mer counts and
GC count on the 200-base core. Their fixed denominators are 198 and 200 in the
v3 model. The code recomputes those integer counts for every source DNA core,
compares every integer and normalized binary64 feature with v3's retained
snapshot, and groups by the **full integer tuples**, without activity labels.
It verifies that no two different tuples acquire the same SHA-256 fingerprint;
fingerprints in the compact assignment ledger are identifiers, not the basis
of equality. Activity labels are joined only after all DNA feature classes have
been constructed.

The 55,276 records contain 55,276 distinct exact core sequences and 55,275
distinct whole-core reverse-complement identities. The feature counts form
55,275 classes: 55,274 singletons and one class of size two. The sole members
are:

| Source ID | hg38 interval, zero-based half-open | Declared/reference-match strand |
| --- | --- | --- |
| `ENSG00000144566` | `chr3:19946928–19947128` | `+` / `+` |
| `ENSG00000163576` | `chr3:19946928–19947128` | `-` / `-` |

Their **200-base cores** are distinct and exactly reverse complementary. The
independent source check also finds that both 230-base S7 oligos have the same
15-base left adapter (`AGGACCGGATCAACT`) and the same 15-base right adapter
(`CATTGCGTGAACCGA`); the complete adapter-containing oligos are not whole
reverse complements. Both records belong to the same previously frozen
dependency component and test fold 4, using zero-based fold numbering (the
fifth of five folds). They are the same genomic interval, not two independently
sampled loci. Source strand annotations and the core comparison do not by
themselves establish the cause of any measured difference.

## Observed-label floor

All numbers below use the first ID minus the second ID. Displayed values are
rounded; retained exact fractions are computed from the actual v3 JSON-decoded
binary64 labels, not from an assumption that their printed decimal spellings are
exact real numbers. The source workbook's XML decimal lexemes may differ from
the JSON spellings even where both decode to the identical binary64 value.

| Cell | First label | Second label | Difference | Pair MSE floor | Whole-panel MSE floor |
| --- | ---: | ---: | ---: | ---: | ---: |
| HepG2 | 1.838 | 1.290 | 0.548 | 0.075076 | 0.00000271640 |
| K562 | 2.204 | 1.876 | 0.328 | 0.026896 | 0.000000973153 |
| WTC11 | 3.142 | 2.868 | 0.274 | 0.018769 | 0.000000679101 |

For a cell with observed labels `y1` and `y2`, any deterministic predictor
using the same complete input for both cores outputs a common `t`. The mean of
their two squared residuals is at least `(y1 - y2)^2 / 4`. The v3 grouping
places both records in the same test fold, so the same per-fold predictor is
used for both when a v3 model uses these features and that cell context. The
whole-panel cell floor weights each of its 55,276 records equally; pooled v3
weighting then averages the three cell floors. The numerator is the sum of
within-class squared deviations, so singleton classes contribute zero. This is
also the empirical optimum of an unrestricted lookup by exact feature vector
and cell on these labels, not a lower bound on every other kind of model.

Only one feature class contributes to this exact-collision floor. The small
whole-panel value does not explain the published model's total error, rule out
other representation limits, or establish the accuracy of predictions on new
sequences. A predictor supplied with orientation, longer-range features, source
ID or other information need not return the same value. The observed differences
can include technical noise, design and assay effects; they do not by themselves
identify a strand-dependent mechanism, an endogenous effect or independent
biological replication.

## Evidence and replay

[`results/assignments.jsonl.gz`](results/assignments.jsonl.gz) retains each
record's exact-feature fingerprint, class size, component and test fold; every
class membership is reconstructible. [`results/collisions.json`](results/collisions.json)
contains the full colliding class, its 33 integer counts, source coordinates,
strand metadata, pair relations, binary64 labels and exact rational cellwise
differences and sums of squares. [`results/summary.json`](results/summary.json)
contains all panel counts and exact rational floors; [`results/manifest.json`](results/manifest.json)
binds the implementation and all three outputs. Rational numerators and
denominators are decimal strings so JSON readers with binary64-only numbers do
not silently round the exact integers.

The implementation verifies the hashes of the panel, v3 feature and activity
payloads, grouping, v3 publication and input receipts, source manifest, frozen
v3 and grouping protocols, and the original source workbook. It checks the
published activity binding and the earlier 258-category independent numerical
reproduction as recorded in the v3 publication. It hashes the raw workbook but
does not reparse its values; it relies on those previously retained source
bindings. The panel JSON lines also contain older ESS values, but this audit
selects and interprets only their DNA and source-metadata fields. The separate
[independent check](../context-activity-feature-collisions-independent-v1/README.md)
reconstructs the result independently and checks the colliding source values.
Hash equality establishes identity to these retained artifacts, not publisher
authentication or empirical truth.

During the initial feasibility inspection, a 25-line preview of the panel file
also printed some previously explored ESS values embedded in those records.
They were not used to identify feature classes, select the pair, or calculate
the activity-label floors; the reproducible program does not interpret them.

From the repository root, select a fresh output directory:

```sh
python3 -m unittest discover -s research/context-activity-feature-collisions-v1 -p 'test_*.py'
python3 research/context-activity-feature-collisions-v1/audit.py --out /tmp/v3-feature-collision-replay
cmp research/context-activity-feature-collisions-v1/results/assignments.jsonl.gz /tmp/v3-feature-collision-replay/assignments.jsonl.gz
cmp research/context-activity-feature-collisions-v1/results/collisions.json /tmp/v3-feature-collision-replay/collisions.json
cmp research/context-activity-feature-collisions-v1/results/summary.json /tmp/v3-feature-collision-replay/summary.json
cmp research/context-activity-feature-collisions-v1/results/manifest.json /tmp/v3-feature-collision-replay/manifest.json
```

The destination must not exist. Tests cover exact feature arithmetic,
reverse-complement invariance, a forced fingerprint collision, binary64 versus
decimal rational arithmetic, weighting and group/fold integrity. This replay
uses the retained raw workbook for hash verification, not an hg38 download or a
refit.
