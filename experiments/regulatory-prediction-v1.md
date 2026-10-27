# Regulatory prediction v1: frozen exploratory specification

The executable authority is
[regulatory-prediction-v1.json](regulatory-prediction-v1.json), SHA-256
`c8e97fabef2afe8ccaa5a01a7ef7bb952f9c9f15857e5c540af112a66e78cb58`.
The final specification was fixed before new score acquisition or real-data
fitting. An initial draft was corrected before acquisition to retain incomplete
positions and use mutation-chemistry-stratified rotations. Further changes require
a new protocol version and an explicit explanation; failed analyses remain in
the record.

This is retrospective exploratory research on public data. The HBB assay and
its historical predictions have already been explored. The newly selected
elements were chosen from metadata before their outcomes entered this workflow,
but they are not an independently preregistered prospective cohort.

## Goal and cohort

Test whether fixed local sequence and prior-motif features add reproducible
information about measured reporter-expression effects beyond mutation identity
and coordinate baselines. Seek candidate explanations that survive controls,
source investigation and independent numerical implementation, while reporting
negative results.

All sources come from Kircher et al. 2019
([DOI:10.1038/s41467-019-11526-w](https://doi.org/10.1038/s41467-019-11526-w)).
Eligibility is human promoter DNA, the same log2 additive regression endpoint,
verified reference identity and CC0 score-set rights. Selection uses no scores.

| Role | Elements / records |
| --- | --- |
| Primary | F9 15-a-1, FOXE1 16-a-1, GP1BB 17-a-1, HBG1 19-a-1, HNF4A 20-a-1, LDLR 23-a-1, MSMB 24-a-1, PKLR 27-a-1, TERT 31-a-1 |
| Validation only | LDLR 23-a-2; PKLR 27-b-1; TERT 31-b-1, 31-c-1, 31-d-1 |
| Previously exposed reference | HBB 18-a-1 |

The table abbreviates full MaveDB URNs, which the JSON lists exactly. For a locus
with several eligible records, the lexicographically first full URN is primary.
Other selected contexts/libraries stay out of training and primary held-element
counts. Cell lines and experimental relationships remain explicit metadata.

Different loci from one study are not independent laboratories or a random
sample of human promoters. Cell context and element identity are confounded in
parts of this cohort. A shared library or a technical barcode does not create an
independent experimental locus.

## Target accounting

Use only unambiguous reference-validated SNVs with finite effects, at least ten
unique tags, and finite positive DNA/RNA counts. No coefficient or p-value cutoff
is used. Missing/invalid/low-tag rows receive explicit exclusions before metrics.
Keep every eligible allele even when another alternate at its position is
missing; report positions with one, two or three eligible alleles.

All alleles at a position share folds, exclusion guards and resampling blocks.
Equality/deletion records remain raw source evidence outside prediction fitting.
No source p-value, read count, tag count, clinical classification, cell identity
or test-set outcome normalization is a sequence feature.

## Fixed features and models

Coordinates refer to the supplied assay target orientation. Only windows wholly
inside that target are used; edge truncation is fixed and disclosed.

| Model | Features beyond an unpenalized intercept |
| --- | --- |
| mean | None |
| chemistry | Twelve directed reference-to-alternate identities |
| coordinate | Chemistry plus normalized target position and its square |
| delta3mer | Chemistry plus 64 signed local 3-mer count changes |
| motif | Chemistry plus four signed motif-score changes |
| full | Chemistry, position, 3-mer changes and motif changes |

The coordinate is `x=2*(position-1)/(length-1)-1`. A delta-3-mer feature counts
alternate minus reference occurrences over the identical fully contained 3-mer
windows that overlap the changed base. Its maximum support radius is two.

Four outcome-independent JASPAR priors are fixed: GATA1 `MA0035.4`, KLF1
`MA0493.3`, NFYA `MA0060.3` and TBP `MA0108.3`. Add a total pseudocount of one per
column, uniformly distributed across bases, and compute log2 odds against uniform
background. Scan every fully contained motif window overlapping the variant on
both strands. The feature is maximum alternate score minus maximum reference
score over the same candidate-window set. Reference maxima are descriptive
outputs, not fitted covariates. No binding threshold is selected.

The widest frozen motif is 11 bp: all affected windows together have radius ten
and maximum span 21 bp. TBP lacks source taxonomy metadata; its interpretation
must retain that limitation rather than inventing a human-derived assay.

Every non-mean model uses fixed ridge alpha **10**, minimizing the sum of squared
errors plus alpha times the squared coefficient norm, without division by sample
size. Training-only feature means and population SDs are used; SD at most
`1e-12` becomes scale one. Center training outcomes by their own row mean, leave
the intercept unpenalized and retain raw log2 outcome units. No per-assay test
intercept, nested hyperparameter search or feature selection is allowed.

## Evaluation that answers distinct questions

**Primary:** nine leave-element-out folds. Train on eight primary loci and test
the remaining locus. HBB and validation-only records never enter those fits.
Average allele errors within a position, positions within an assay, then assays
equally. Report each element and context as well as the macro result.

**Exploratory HBB:** five contiguous-position folds with boundaries
`floor(k*length/5)`. Remove training centers whose distance from any test center
is at most **20 bp**. Radius-ten feature supports then cannot overlap between
training and test centers. Source dependence and regional biology can remain.

**Transfer diagnostic:** fit on all eligible HBB records only and apply the frozen
transform to all nine primary elements without test normalization.

**Explicitly leaky diagnostic:** HBB-only row-random five-fold assignment from the
specified SHA-256 rule, with no grouping guard. Report shared positions and
overlapping sequence support. A gap from guarded evaluation also reflects a
different split difficulty; it is not a causal estimate of leakage alone.

Primary error is MSE in squared log2-effect units. Secondary outputs are MAE,
row-level Spearman, position-mean Spearman and per-assay relative MSE against the
training-mean baseline. Spearman uses average ranks for exact ties and is
unavailable, rather than zero, for constant vectors.

## Prespecified comparisons and controls

There are exactly three candidate-minus-baseline comparisons:

1. delta3mer versus chemistry.
2. motif versus chemistry.
3. full versus coordinate.

Improvement is baseline MSE minus candidate MSE; larger is better.

Use 99 spatial rotation controls in both the guarded-HBB and primary
leave-element-out families. Within each assay and directed mutation class, order
eligible rows by position and circularly rotate their effects by a nonzero
offset. Singleton strata remain unchanged and are disclosed. All features,
eligible IDs and folds remain fixed. Recompute the rotated training mean and
right-hand side, reusing fixed feature factorizations for efficiency.

This preserves mutation-class distributions and cyclic positional ordering
within a class. It breaks some cross-class/within-position correlation and does
not preserve absolute genomic gaps or all library dependence. Spatial stationarity
and exchangeability are not established for promoters. Therefore empirical
rotation tail fractions are **diagnostic controls, not calibrated
biological-significance p-values**.

Report the plus-one tail fraction, Holm adjustment across the three comparisons
within each family, and conservative Bonferroni adjustment over all six planned
tails. Do not select the better family, model or null after viewing results.

Two thousand common-draw resamples summarize conditional stability of frozen
out-of-fold losses:

- HBB: nonoverlapping 21-bp physical position blocks, including the shorter last
  block; draw the same number of nonempty blocks with replacement.
- Primary cohort: nine element indices drawn with replacement, retaining equal
  element weighting.

All three comparisons use the same draws within a resampling family. The JSON
fixes SplitMix64 arithmetic, seeds and Type-7 percentile quantiles. These are
descriptive intervals conditional on fixed predictions, not simultaneous
confidence guarantees or uncertainty over fitting/model selection. Overlapping
out-of-fold fits and one shared source study limit inference.

## Validation-only context and reliability analysis

Match exact genomic SNVs across the prespecified LDLR library repeat, PKLR
timepoints and TERT cell/perturbation contexts. Report overlap, exclusions,
correlation, MSE/MAE and signed mean differences without fitting on these records.

For a variant observed in several contexts, any predictor assigning one
context-blind value has observed mean-square error at least

`sum_variant,context (observed_effect - variant_mean_effect)^2 / observations`.

This is a finite algebraic lower bound on these measurements. It does not
attribute variation causally to cell biology or separate biological variation
from noise, normalization and library differences. Compare technical/library
differences with context differences while retaining their source genealogy.

## Reproduction and claims

Independent Rust and Python implementations must agree on included IDs,
exclusions, feature vocabulary, folds/guards, training scaling and predictions.
Prediction absolute tolerance is `1e-7` and metric tolerance `1e-8`. Investigate
failures; do not loosen tolerances silently.

Each new output directory retains the protocol/input hashes, full predictions,
coefficients/scaling, diagnostics, controls and failures. A pattern is at most a
candidate mechanism when it survives relevant controls and source scrutiny.
Check prior literature before describing novelty; rediscoveries and curation
corrections are useful but must be named accurately. Independent perturbation or
laboratory data are still required for a mechanistic claim.
