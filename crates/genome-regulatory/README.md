# Frozen regulatory-prediction engine

This crate executes `experiments/regulatory-prediction-v1.json`, pinned by its
final pre-acquisition SHA-256. It accepts reporter-effect records, not clinical
classifications. The study is retrospective and HBB outcomes were already
explored; it does not claim pristine HBB validation or a complete genome model.

```rust,ignore
let protocol = genome_regulatory::Protocol::open(protocol_path)?;
let motifs = genome_regulatory::motif::MotifSet::open(motif_directory)?;
let study = genome_regulatory::run(&records, &protocol, &motifs)?;
```

For a complete standalone reproduction, run from the repository root:

```bash
source scripts/env.sh
cargo run --release --locked -p genome-regulatory --example frozen_study -- \
  /tmp/a-new-regulatory-study.json
```

Output files must not already exist. The example loads the frozen regulatory
cohort, which already includes HBB; do not append HBB a second time.

## Data and eligibility

`InputRecord` holds an accession, assay ID, one-based target position,
reference/alternate base, target sequence, cell-context description and optional
score/tag/DNA/RNA values. Assembly and target-relative HGVS reconciliation are
the responsibility of the source adapter; the engine verifies target reference
bases and unique accessions/SNVs again.

`score_raw` also retains the original CSV decimal, whose parsed value must agree
with the supplied numeric outcome. Observed position means and their rank ties
use exact rational arithmetic on these literals. The small `num-bigint` and
`num-rational` dependencies prevent an order of floating-point summation from
turning equal source means into different ranks. Predicted values are never
quantized or rounded to manufacture ties.

Quality requires finite scores, at least ten integral unique tags, and positive
finite integral DNA/RNA counts. Missing values are explicit exclusions; a direct
`Some(NaN)` or infinity is rejected. Each passing allele remains eligible when a
sibling allele is missing. Included IDs, exclusion reasons and one/two/three-allele
position counts are all saved.

Primary training uses only the nine protocol-designated primary loci.
Validation-only records and HBB never enter leave-element-out training. HBB is
used only in its declared exploratory within-element/transport analyses.
Assay/cell IDs and source counts/p-values are never learned covariates.

## Fixed feature definitions

| Model | Features |
| --- | --- |
| `mean` | Training outcome mean only |
| `chemistry` | Twelve directed substitution identities |
| `coordinate` | Chemistry plus normalized target position and its square |
| `delta3mer` | Chemistry plus 64 alternate-minus-reference local 3-mer counts |
| `motif` | Chemistry plus four fixed JASPAR signed log-odds deltas |
| `full` | All of the above, 82 features |

JASPAR scanning uses complete windows overlapping the changed base, on both
strands, with identical candidate windows for reference and alternate.
Probabilities are `(count + 0.25)/(column_total + 1)` against uniform 0.25
background. Reference/alternate maxima and their locations are retained for
inspection; only four signed deltas enter a predictor. Similarity scores are not
binding energies, occupancy measurements or expression coefficients.

Every model uses fixed ridge alpha 10, an unpenalized intercept, and training-only
population-SD feature standardization. SD at most `1e-12` uses scale one. The
objective is summed squared error plus alpha times the squared coefficient norm.
No model/window/alpha selection is performed.

## Evaluation and controls

- HBB: five contiguous position folds, boundaries `floor(k*length/5)`, with
  training centers more than 20 bases from every test center. The maximum
  feature-support radius is ten, so these sequence supports are disjoint.
- Nine primary elements: leave one locus out, training on the other eight.
- Secondary transport: fit HBB and evaluate each primary element unchanged.
- Leaky diagnostic: five deterministic SHA-256 row folds on HBB. Shared
  positions and feature-support overlap are explicitly reported.

MSE/MAE first average retained alleles within a position, then positions within
an assay, then assays equally. Row and position-mean correlations are secondary.
Exact constant vectors yield undefined correlation; constant prediction means
remain exact when position groups contain different allele counts. Negative
out-of-sample R-squared is preserved.

The result retains all targets, predictions, fold/guard IDs, fitted intercepts,
coefficients and training means/scales. It also retains 2,000 paired physical-block
and element bootstrap draws, plus 99 chemistry-stratified spatial-rotation
controls. One fixed SplitMix64 stream rotates scores within assay and directed
mutation class. Every null fit recomputes its outcome mean and right-hand side;
only outcome-independent X preprocessing and Cholesky factorizations are cached.

Bootstrap intervals condition on fixed OOF losses. Rotation tails describe an
artificial spatial null, not calibrated biological significance. Holm adjustments
within three-comparison families and conservative six-tail Bonferroni values are
both reported. Shared studies, cell-context differences and prior exposure remain
limitations; thresholds never turn a pattern into a discovery.

Independent replay reproduced all 131,112 predictions within the fixed numerical
tolerance, as well as folds, fitted parameters, losses, bootstraps and nulls.
Twenty-six secondary exact Spearman values remain implementation-sensitive:
floating-point ridge solutions can split mathematically equal effects by a few
ulps. The replication report retains these mismatches. No epsilon tie threshold,
prediction rounding or tolerance relaxation was introduced to conceal them.
