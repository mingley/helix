# Independent numerical replication

`scripts/replicate_regulatory_study.py` reconstructs the frozen study from raw
MaveDB CSVs, reference strings and JASPAR PFMs using NumPy/SciPy. It does not call
Rust feature/fold/fitting code or use Rust predictions to fit models.

This verifies software. It is not independent laboratory replication or
validation of the entire human genome.

## Frozen authority and environment

The accepted protocol is `experiments/regulatory-prediction-v1.json`, SHA-256:

    c8e97fabef2afe8ccaa5a01a7ef7bb952f9c9f15857e5c540af112a66e78cb58

Real outcome parsing requires this frozen protocol, verified input commitments,
and matching cohort acquisition gates/roles. HBB was previously explored; the
nine other primary promoters were acquired under the frozen plan. They share
one study and source pipeline.

The original run used CPython 3.12.14, NumPy 2.3.5 and SciPy 1.17.0. Each run
records actual versions and retains an exact implementation snapshot. Use the
repository's pinned Python requirements. One BLAS thread avoids competing with
concurrent work; it is a numerical convention, not a model parameter.

## Run and compare

From the repository root, choose new output directories:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  python3 scripts/replicate_regulatory_study.py \
  --out /tmp/independent-regulatory-new

OPENBLAS_NUM_THREADS=1 python3 scripts/replicate_regulatory_study.py \
  --compare-rust /path/to/study.json.gz \
  --independent-run /tmp/independent-regulatory-new \
  --out /tmp/independent-comparison-new
```

The first command executes all six models, four evaluation schemes, 99 cached
spatial-null refits, 2,000 block/element bootstrap draws, and fixed context pairs.
No parameter is selected using outcomes.

Comparison accepts plain or gzip JSON and records both compressed and decoded
commitments. Existing directories are refused. Failed checks are saved and
exit nonzero, including failures in secondary statistics.

## Independent calculations

- Parse original source rows using fixed score/tag/read-count quality rules;
  retain eligible alleles at incomplete positions.
- Reconstruct mutation identities, coordinates, signed overlapping 3-mer counts
  and dual-strand motif deltas directly from DNA/PFM inputs.
- Construct physical folds/guards from positions and raw source IDs.
- Use training-only population mean/SD scaling, unit scale below the frozen
  tolerance, an unpenalized intercept and the fixed summed-error ridge penalty.
- Cache SciPy Cholesky factors, then recompute every rotated training outcome
  mean and right-hand side.
- Independently implement SplitMix64 and shared bootstrap/rotation draw tables.
- Weight errors equally across alleles within position, positions within assay,
  and assays within the primary cohort.

The comparator checks every prediction/fold, feature vocabulary/dimension,
training preprocessing, coefficient/intercept, per-assay loss, bootstrap draw
and statistic, and null contrast. Frozen tolerances remain 1e-7 absolute for
predictions and 1e-8 for metrics.

## Rank hazards and retained attempts

Initial attempts exposed two numerical artifacts. Averaging a constant predictor
over groups with one, two or three alleles can produce one-ULP differences.
BLAS batch alignment can also split predictions of identical feature vectors.
Exact-tie ranks can turn those artifacts into apparent correlations.

Corrections preserve mathematical definitions:

1. Evaluate identical feature inputs once per fit, then broadcast.
2. Rank observed position means using exact rational arithmetic from raw
   decimal score strings.
3. Preserve exact constant groups/vectors and report undefined correlations null.

There is no score quantization, rank epsilon, feature change or model selection.
Distinct inputs can still yield theoretically equal predictions whose independent
floating-point fits differ near machine precision. Remaining exact-rank
mismatches stay failed secondary checks. Successful primary prediction/loss
checks are reported separately and do not establish complete reproduction.

Failed and corrected attempts remain separate bundles. Reports retain
discrepancies, exact commitments and maximum absolute errors.

The [published evidence snapshot](../results/regulatory-v1/README.md) contains
the full final independent bundle, failed comparison and compressed historical
attempts. Earlier absolute `/tmp` paths in the original receipts describe where
those attempts ran; the snapshot's `history/attempts.json` maps them to portable
archives with original/archive hashes. Historical implementation snapshots were
not captured for every early attempt; the index discloses that limit.

## Synthetic checks

```bash
OPENBLAS_NUM_THREADS=1 python3 -m unittest discover \
  -s research/independent-regulatory-replication -p test_replication.py -v
```

Cases cover a closed-form ridge result, training-only scaling, unpenalized
intercept, sparse positions, guard zones, signed k-mers, hand-computed motif
log-odds, reverse-complement invariance, PRNG fixtures, rotation order, cached
refit means, position weighting and numerical rank regressions. No network
or real assay fitting is needed.

Motif similarity is not binding or expression. HBB's physical insert is in
native-minus orientation relative to the coordinate target. Dual-strand motif
deltas are reversal-invariant, while reporter transcript/TSS identity remains
unmeasured. No RNA/uAUG mechanism is inferred.
