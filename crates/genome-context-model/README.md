# Matched-context sequence models

This crate implements five small regression models for the same assayed DNA
sequences measured in multiple observed experimental contexts. It supplies
sequence features, complete-panel validation, grouped folds, training-only
ridge fits, held-out predictions, paired context contrasts, and reproducible
component bootstrap summaries. Source selection, endpoint construction, reuse
rights, grouping, fold assignment, fixed hyperparameters, and scientific
interpretation belong to the separately frozen study protocol.

The implementation was prepared and tested using synthetic outcomes before
accessing the proposed next study's outcomes. Creating this crate does not
establish that a candidate source contains the required observed endpoint.
Predicted source scores, expression-specificity summaries, and clinical labels
must not be passed off as measured reporter activities.

## Public contracts

```rust,ignore
use genome_context_model::{
    evaluate, paired_group_bootstrap, GroupSplit, LossWeighting,
    MatchedPanel, MeasuredOutcome, ModelComparison, ModelKind,
    OutcomeMatrix, SequenceRecord,
};

// Sequence-only metadata; grouping is computed independently of outcomes.
let panel = MatchedPanel::new(
    vec![/* SequenceRecord { id, group, sequence }, ... */],
    vec!["context-a".into(), "context-b".into()],
)?;
let folds = GroupSplit::partition(&panel, &fold_test_groups)?;
// For a prespecified physical-neighborhood guard, construct each split with
// GroupSplit::with_guard(&panel, &test_groups, &excluded_training_sequence_ids).

// Explicit keys avoid ambiguous matrix/column alignment. Every sequence must
// have exactly one finite outcome in every context, including held-out rows.
let matrix = OutcomeMatrix::from_records(
    &panel,
    &outcome_rows, // Vec<MeasuredOutcome { sequence_id, context, value }>
)?;
let evaluation = evaluate(&panel, &folds, &matrix, 10.0)?;

let comparisons = [
    ModelComparison { baseline: ModelKind::SharedDna,
                      candidate: ModelKind::AdditiveContextDna },
    ModelComparison { baseline: ModelKind::ContextMean,
                      candidate: ModelKind::AdditiveContextDna },
    ModelComparison { baseline: ModelKind::AdditiveContextDna,
                      candidate: ModelKind::IndependentContextDna },
];
let bootstrap = paired_group_bootstrap(
    &evaluation.per_group_errors, &comparisons,
    2_000, 0x2026102702, [0.025, 0.975], LossWeighting::EqualSequences,
)?;
let bonferroni_three_ranges = bootstrap.intervals_at([0.05/6.0, 1.0-0.05/6.0])?;
```

`MatchedPanel` sorts unique sequence IDs and context names and rejects duplicate
IDs, empty/padded labels, non-ACGT DNA, or sequences shorter than three bases.
Exact and reverse-complement-identical DNA may have multiple source IDs, but
those IDs must share a declared leakage group. This check alone does not
establish family, homology, genomic-distance, or source-parent independence.
Those rules must be validated from independently pinned metadata.

`OutcomeMatrix` rejects unknown keys, duplicate measurements, incomplete pairs,
and nonfinite values. It does not impute, transform, normalize, average technical
replicates, interpret missing-value sentinels, or choose an endpoint. A source
loader must preserve excluded rows and explain any prior aggregation and
source-level normalization. These are measurements of the supplied assay
endpoint, not automatically endogenous regulatory activity.

`GroupSplit::new` holds out complete groups in every context.
`GroupSplit::with_guard` additionally excludes a metadata-derived list of
training sequence IDs in every context; test groups remain intact. A protocol
requiring a whole-component guard must supply all IDs from each excluded
component. Splits preserve train/test/guard IDs, groups, and contexts.
`evaluate` requires every sequence to occur in exactly one test fold. Split and
outcome objects are bound to the panel's metadata digest, so they cannot be
reused against a different panel accidentally.

## Features and models

There are 33 outcome-independent features: 32 reverse-complement-canonical
trinucleotide frequencies in lexicographic A/C/G/T order, followed by GC
fraction. Each overlapping complete 3-mer contributes one integer count to
`min(kmer, reverse_complement(kmer))`; divide each count once by `length-2`.
GC is `(G+C)/length`. Feature names are `rc3mer_AAA` through `rc3mer_TCA`, then
`gc_fraction`. The full ordered list is returned by `feature_names()` and
stored in every fit. Sequence length and raw counts are available for audits;
length is not an extra predictor. Sequence features are invariant to reverse
complementation and contain no context, outcome, count-quality or clinical
information.

Each fold prepares one 33-feature design using **unique training sequences**.
Population SD uses denominator `n`; SD at most `1e-12` is replaced by one. The
same training means, scales, and Cholesky factorization are reused for shared
and context-specific fits. Test DNA does not influence this transformation.
The existing `genome-regulatory::math::RidgeDesign` supplies the deterministic
small ridge solve. No context-by-DNA dense expanded matrix is allocated.

| Model | Training target and prediction |
| --- | --- |
| G, `global_mean` | Grand training mean, shared by all contexts |
| C, `context_mean` | Training mean for the known context |
| S, `shared_dna` | Ridge on one context-mean target per training sequence |
| A, `additive_context_dna` | The exact S slopes, with free training context intercepts |
| I, `independent_context_dna` | One ridge fit per context on the same training DNA/scaler |

The objective for S is summed squared error over `n` sequence-mean targets plus
`alpha * ||beta||²`. Its equivalent `n*C`-measurement pooled objective has
penalty `C*alpha`, not `alpha`. Each I fit uses `n` training targets and the same
alpha. Intercepts are unpenalized and outcomes remain on their supplied scale.
The library requires an explicit positive alpha; an experiment must fix it
before evaluating outcomes.

A predicts `training_context_mean + shared_DNA_dot_product` directly, avoiding
cancellation through the grand mean. On a complete panel, A and S have the
same slopes. Linearity also gives `mean_context(I_beta) = S_beta` and mean-context
I predictions equal S predictions, up to floating-point reduction errors.
The output reports these numerical checks. The synthetic suite also verifies
that changing held-out outcomes cannot alter fits or predictions.

## Losses, contrasts and retained artifacts

`ContextEvaluation` serializes all fold memberships and fitted coefficients,
means/scales, training-only outcome digests, every paired held-out prediction,
per-sequence/context/group MSE and MAE, and both sequence-weighted and
macro-group losses. All five models have identical target measurements.
The caller must identify its primary estimand. Equal sequence weighting and
equal component weighting differ when component sizes differ.

All context pairs have right-minus-left observed and predicted differences,
plus per-component errors. On a complete panel, the saved checks verify

```text
MSE(A) - MSE(I)
  = (C-1)/(2*C) * mean_pair[contrastMSE(A) - contrastMSE(I)].
```

This is an algebraic check of predictive interactions, not an additional
independent experiment. Known-context offsets do not establish context-specific
DNA effects; the A versus I comparison addresses that narrower model question.

`paired_group_bootstrap` samples the same number of components with replacement
using one `SplitMix64` stream. Group order is lexicographic; all comparisons,
contexts, and MSE/MAE summaries share the exact draw table. `EqualSequences`
weights a sampled component by its sequence count and divides pooled sampled
loss sums by pooled sampled sequence counts. `EqualGroups` instead averages
component means and is a distinct sensitivity estimand. All indices and sampled
statistics are returned. Additional prespecified percentile levels reuse these
same statistics via `intervals_at`, without repeating draws or fitting.
All five models also retain observed and sampled MSE/MAE arrays, with explicit
`model_order`, plus per-model percentile summaries. They reuse the same component
indices and weights as the comparisons and consume no additional random draws.
`intervals_at` returns both comparison and model ranges from those saved values.
Quantiles use Type-7 linear interpolation. Resource guards reject more than
100000 replicates or more than 50 million retained sample indices.

Panel and outcome digests commit the canonical, typed inputs; the caller must
also pin raw source files, source/endpoint provenance, protocol, grouping,
software revision, and artifact hashes. Conditional bootstrap ranges omit
model-refitting uncertainty and do not create independent studies or imply
genome-wide generalization. Composition models do not identify transcription
factor binding, causal molecular mechanisms, clinical outcomes, or new contexts.

## Verification

```sh
source scripts/env.sh
cargo test -p genome-context-model --offline
cargo clippy -p genome-context-model --all-targets --offline -- -D warnings
```

Tests cover a closed-form ridge value, exact overlapping 3-mer counts and
reverse-complement invariance, complete-pair/duplicate/nonfinite rejection,
group and guard isolation, held-out-outcome perturbation, common training-only
scaling, constant outcomes, S/A/I linearity, the three-context contrast identity,
unequal component weighting, and shared deterministic bootstrap draws.
