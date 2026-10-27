# Draft: context-aware prediction for held-out regulatory sequences

**Status: design draft, not an executable frozen protocol.** This document was
written before acquiring outcomes for the proposed next study. It authorizes no
outcome download or fitting. Source accessions, permission to reuse score data,
endpoint columns, missingness rules and final grouping feasibility must be
verified and frozen separately before acquisition. The previous Kircher/HBB
results remain exploratory and do not become a new holdout.

## Question and scope

Can observed experimental context improve prediction of a new regulatory
sequence's measured activity, after accounting for both the context's average
activity and a shared DNA model? Does the relationship between DNA composition
and activity vary across contexts beyond an additive context offset?

The target is **new sequence families in already observed assay contexts**.
Cells, vector designs, integration status, insert location, time and perturbation
are distinct metadata fields. A difference between two complete assay contexts
cannot by itself be attributed to one biological factor. This study does not
test an unseen cell type, a new source study, variant pathogenicity or an
unselected human genome. Whole-element activity is a different endpoint from
the previous single-nucleotide effect coefficients.

The source-search team has nominated Agarwal et al. (2025),
DOI `10.1038/s41586-024-08430-9`, with a reportedly common approximately 60000
cCRE panel in HepG2, K562 and WTC11. This nomination is metadata-level only in
this draft; exact files and reuse rights are pending verification. Prior
activity-based selection and additional variable-promoter selection mean that
such a panel is an experimentally selected population. Results would concern
that panel, even if the current score table remains unexamined until freezing.
Klein/Inoue multi-design panels are alternatives only if they meet the same
endpoint, paired-sequence, mapping and licensing gates. A source will not be
selected because of its observed predictability or apparent context effects.

For the nominated source, the proposed endpoint is the mean of published
per-replicate median-normalized log2 RNA/DNA activities, without an additional
across-panel z-score. The source paper already studies sequence/context
prediction and uses a full-panel per-cell z-scoring/averaging construction for
another analysis. This small fixed-model, grouped evaluation would therefore
test a distinct explicit endpoint and control structure; the general idea that
regulatory activity depends on context is not a novel finding. The source
team must verify the exact replicate fields and formula before freezing them.

## Source gate before outcomes

Freeze a source ledger with publication, persistent accession, source genealogy,
raw URLs and hashes, score-data license scope, reference assembly, insert
sequences, sequence-to-parent mapping, coordinate convention, context factors,
replicate design, source endpoint and its units. Metadata-only licensing is
insufficient permission to redistribute scores. Distinguish independent source
studies from reuploads, shared libraries, technical replicates and reused
sequences. A second cell context within the same study is not a second study.

The primary panel must contain the same actual inserted DNA for each sequence
identity in every chosen context. A shared genomic parent with different-length
inserts is insufficient for the balanced-panel algebra below. Verify this from
sequence assets, rather than treating a common label as proof. Context-specific
insert sequences require a different prospective protocol.

Choose the documented activity summary and transform from assay methods before
reading its values. Prefer a published replicate-aggregated endpoint with clear
normalization and finite missing-value conventions. Replicates, tags, barcodes,
RNA/DNA molecules and repeated source rows are not independent sequences.
Predeclare aggregation or retain the published aggregate; do not choose an
aggregation after seeing its performance. Do not filter by measured activity,
activity ranks, source p-values or observed replicate agreement. Any necessary
source count/QC rule must be explicit and its possible activity-dependent
missingness disclosed.

Use complete sequence-by-context rows for the primary analysis. Preserve all
incomplete records and a context-specific exclusion table for a separate
coverage report. No imputation from other contexts and no favorable subset
selected after fitting. Complete-case inference concerns complete cases, not
every designed insert. Freeze the required number of matched contexts and
minimum sample/group counts after metadata inspection but before score access.

Source-provided normalization may already have used the entire experiment.
Record this provenance: training-only processing in our pipeline cannot undo
global source normalization. Additional centering, scaling or transformation
learned by our models must use training data only. No test-context mean or test
sequence activity may enter a fit.

## Outcome-independent sequence groups and folds

The proposed default is five deterministic outer folds of indivisible sequence
groups. All contexts of a sequence use its one fold. A random row split would
place the same insert in training and test and is not an admissible primary
evaluation.

Build an undirected union graph using metadata and DNA only. Join sequences
sharing any of the following:

1. A documented source parent/family identity.
2. Exact sequence or reverse-complement identity.
3. Positive-length genomic overlap on the same verified assembly and chromosome.
4. A touched 1000000-bp genomic bin, anchored at coordinate zero. A half-open
   interval `[start,end)` touches bins `start//1000000` through
   `(end-1)//1000000`; an interval crossing a boundary joins both bins.
5. At least 90% global, ungapped identity for equal-length inserts in either
   orientation: minimum forward/reverse-complement Hamming distance at most
   `floor(length/10)`.

Connected components are indivisible, including transitive chains. The final
metadata freeze must report component sizes, genomic coverage and sources of
each union before assigning folds. A proposed feasibility gate is at least 25
components, at least five per test fold and no component larger than 20% of
sequences. If a metadata-only gate fails, stop and amend the draft before any
outcome acquisition; do not silently split a component or weaken identity rules
to obtain a favorable result.

The equal-length similarity rule can be exhaustive without an all-pairs scan.
For each length `L`, let `k=floor(L/10)` and partition the full sequence into
`k+1` nonempty contiguous segments using boundaries `floor(j*L/(k+1))`.
Index both orientations by `(L,segment_index,exact_segment)`, deduplicate
candidate sequence pairs, verify their full Hamming distance and union only
verified pairs. At most `k` mismatches leave at least one of these `k+1`
segments unchanged. This is a complete candidate procedure for the stated
ungapped rule; approximate sketches alone do not supply that guarantee. Do not
silently truncate large low-complexity buckets. Exact duplicate collapsing
before candidate generation is permissible if all source IDs remain recorded.
Unequal-length homology, indels and distant shared regulatory motifs are not
ruled out by this rule and remain limitations.

A computational budget may stop this metadata-only step, but must not replace
the exact rule with dropped buckets, a reduced favorable panel or approximate
matches described as exhaustive. Record an abort and fix the implementation or
amend the design before outcome acquisition. Freeze any numeric budget with the
metadata-only implementation.

Sort components by decreasing unique-sequence count, then ascending SHA-256 of
their sorted stable source IDs joined by a newline (IDs must not contain a
newline). Assign each component to the fold with the smallest current sequence
count, breaking ties by fold index 0 through 4. Record membership and the
algorithm version, rather than relying on process hash order.

One-megabase bins alone allow close neighbors on opposite bin boundaries.
Proposed additional guard: for a held-out fold, exclude a training component if
any of its verified intervals is within 10000 bp of any test interval on the
same chromosome. Interval distance is zero for overlap or adjacency, otherwise
the intervening half-open gap; distance at most 10000 is excluded. Exclude the
whole component and every context together. Keep all test sequences. Validate
the final minimum distance and fold training counts before outcomes. The guard
is a fixed spatial robustness choice, not a claim that all correlations end at
10 kb. Its final inclusion and exact feasibility must be frozen with the source
coordinates.

## Fixed DNA features and five models

Proposed accepted model family: 32 reverse-complement-canonical trinucleotide
frequencies plus GC fraction, 33 features in total. Enumerate canonical strings
lexicographically over A/C/G/T, map each overlapping 3-mer to the smaller of it
and its reverse complement, count in integers and divide once by `L-2`. GC is
`(count(G)+count(C))/L`. Require documented insert DNA with only A/C/G/T and
length at least three. Do not include vector flanks, context IDs as DNA,
activity, tag counts, source p-values, replicate outcomes or clinical labels.
This modest composition model cannot resolve motif order, chromatin state or
every context-dependent mechanism. A negative result is limited to this model.

The exact proposed order is the lexicographically sorted canonical 3-mers
followed by GC:

```text
rc3mer_AAA rc3mer_AAC rc3mer_AAG rc3mer_AAT rc3mer_ACA rc3mer_ACC
rc3mer_ACG rc3mer_ACT rc3mer_AGA rc3mer_AGC rc3mer_AGG rc3mer_ATA
rc3mer_ATC rc3mer_ATG rc3mer_CAA rc3mer_CAC rc3mer_CAG rc3mer_CCA
rc3mer_CCC rc3mer_CCG rc3mer_CGA rc3mer_CGC rc3mer_CTA rc3mer_CTC
rc3mer_GAA rc3mer_GAC rc3mer_GCA rc3mer_GCC rc3mer_GGA rc3mer_GTA
rc3mer_TAA rc3mer_TCA gc_fraction
```

For each fold, calculate feature means and population SDs on **unique training
sequences**, once. Set a scale to one when its SD is at most `1e-12`. Apply this
same training transform to all contexts and test sequences. Use fixed ridge
`alpha=10`, a centered response and an unpenalized intercept; no tuning or
feature selection. All models use identical eligible sequences and folds.

Let training sequences be `s=1..n`, observed contexts `c=1..C`, standardized
features `z_s`, context means `m_c=mean_s(y_sc)`, sequence means
`a_s=mean_c(y_sc)` and grand mean `m=mean_c(m_c)`.

| Model | Training and held-out prediction |
| --- | --- |
| G: global mean | `m` |
| C: context means | `m_c` |
| S: shared DNA | Fit `beta=(Z'Z+10I)^-1 Z'(a-m)`; predict `m+z*beta` |
| A: additive DNA + context | Reuse the exact S slope object; predict `m_c+z*beta` |
| I: context-specific DNA | For each context fit `beta_c=(Z'Z+10I)^-1 Z'(y_c-m_c)`; predict `m_c+z*beta_c` |

The shared objective is sum-squared error on `n` sequence-mean targets plus
`10*||beta||^2`. Equivalently, its pooled `n*C`-row objective uses penalty
`C*10*||beta||^2`. Independent context fits use `n` rows and penalty ten each.
Using penalty ten on pooled shared rows instead would change effective
regularization with context count; it is not this design. All models retain
native activity units; there is no learned test-set normalization.

The balanced design supplies useful falsification checks. S and A must have
identical slopes; their predictions differ exactly by the training context
mean minus the grand mean. Linearity of ridge in the response also gives
`mean_c(beta_c)=beta`, and hence mean-context I predictions equal the S/A
mean-context prediction, up to a fixed numerical validation tolerance. These
invariants should pass synthetic tests and the saved fit-artifact checks. If
they fail, investigate implementation or missing balance before interpreting
results.

Compute A predictions by adding the stored shared-DNA dot product directly to
the training context mean. Independent implementations may use different linear
algebra reductions, so exact membership/feature checks and fixed numerical
tolerances for coefficients, predictions and losses remain distinct. Proposed
prediction and loss tolerances are `1e-7` and `1e-8`; do not widen them after
seeing a failed comparison.

I has `C` slope vectors and greater flexibility than A. The held-out comparison
penalizes unsuccessful flexibility, but a positive result supports this
particular context-dependent composition model rather than identifying a
causal molecular interaction. Context identities are observed labels, not
features enabling prediction for an unseen context.

## Primary losses, contrasts and uncertainty

For each out-of-fold sequence, average squared error across the fixed contexts;
then average sequences equally. On a complete panel this also equals an equal
average of the per-context sequence MSEs. Report native squared activity units,
MAE, all context-specific errors, all five model errors and eligible counts.
MSE is primary; MAE is a prespecified robustness endpoint. Do not use exact-rank
correlation as an acceptance criterion. Predictions are cross-fitted once;
held-out outcomes are used only for scoring.

The three prespecified improvements are baseline loss minus candidate loss:

| Contrast | Falsifiable question |
| --- | --- |
| S minus A | Does known context improve shared-DNA prediction on new sequences? |
| C minus A | Does DNA improve beyond context offsets? |
| A minus I | Do context-specific sequence slopes improve beyond additive offsets? |

All three must be retained, favorable or unfavorable. A favorable S-minus-A
contrast alone demonstrates predictable offsets, not context-dependent DNA
effects. The strongest intended claim requires improvements beyond both simple
controls and the interaction comparison. Report absolute improvements and
relative improvements with the denominator named; do not elevate a small
selected difference into genome-wide mechanistic discovery.

Use 2000 common-draw cluster bootstrap resamples of the indivisible genomic/
sequence-family components. Keep every context of each sampled sequence
together. For a draw, sum sequence-level paired loss differences in the sampled
components, with multiplicity, and divide by the number of sampled sequences.
Do not equally average differently sized component means unless changing the
estimand explicitly. Use the same draws for all three contrasts. Report
ordinary 95% conditional ranges and, if adopting simultaneous intervals at
freeze, Bonferroni three-comparison ranges with quantiles `0.05/6` and
`1-0.05/6`; freeze Type-7 interpolation and the PRNG/seed before running.

These are conditional resampling summaries for fixed out-of-fold predictions
on the selected source panel. They omit refitting uncertainty, do not create
independent studies and do not estimate generalization to arbitrary human
genomic sequences. Genomic component bootstrap is more defensible than treating
barcodes or sequence-context rows as independent, but cannot prove residual
component independence. No artificial permutation tail or calibrated p-value
is required by this draft. Any later null control is a separate prespecified
diagnostic, not retroactive evidence of significance.

For every context pair, report MSE/MAE of the observed activity difference and
its predicted difference using A and I. A predicts the training context-mean
difference; I additionally predicts sequence-dependent differences. With the
same complete panel and weights,

`MSE(A)-MSE(I) = (C-1)/(2*C) * mean_pair[contrastMSE(A)-contrastMSE(I)]`.

This follows because mean-context predictions agree and pairwise differences
capture the centered context component. It is an artifact check and a useful
description of interaction, not another independent experiment or an extra
selected positive endpoint. A table of all pairs is required; no best cell pair
is selected after scoring.

## Acceptance, reporting and remaining freeze decisions

Before scoring, freeze raw asset hashes, inclusion/exclusion rules, the actual
insert/context matching audit, grouping/fold/guard algorithms, all five model
definitions, feature order, ridge/scaling conventions, metrics, bootstrap seed
and quantiles, the three contrasts and numeric replication tolerances.
Independently reconstruct grouping and feature vectors and independently fit
the models. Compare saved memberships, training-only transforms, coefficients,
out-of-fold predictions and losses. Deliberately perturb test-only outcomes in
a synthetic test and require every fitted quantity to remain unchanged. Tests
also need RC feature invariance, overlap/boundary grouping, near-identity
pigeonhole completeness on small exhaustive examples, incomplete-panel
rejection and the balanced-model identities above.

Publish a full source-selection and attempt ledger, QC/missingness by context,
group sizes, fold balance/guard exclusions, all model losses, all planned
contrasts, all pairwise context diagnostics and independent-comparison status.
Retain failed attempts and negative findings. Any outcome-informed amendment is
a separately versioned exploratory supplement, with the original result kept.

Remaining source-dependent decisions are explicitly unresolved here: final
source files/license, whether matching actual inserts is possible, endpoint and
published normalization, replicate aggregation, QC/missingness conventions,
sequence/coordinate validity, whether 1-Mb/identity grouping passes feasibility,
whether the 10-kb whole-component guard is retained, minimum primary support,
bootstrap seed and simultaneous-range convention. None should be resolved by
looking for a better outcome or score.
