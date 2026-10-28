# A design for testing the model's predicted zero interaction

This is a **source-agnostic methods note**, not a frozen analysis protocol or a
biological result. The next step is to identify a public assay whose actual DNA,
activity scale and replication satisfy the requirements below, and then freeze
the source-specific cohort and analysis before inspecting its outcomes. These
requirements are distinct from the exact model result already established in
the [locality study](variant-locality-study.md).

## The square and its measurement

For one reference insert and one cell context, define four same-length DNA
constructs: `x00` is the reference, `x10` and `x01` contain one substitution
each, and `x11` contains exactly both substitutions. At each edited position,
the single and double constructs must have the same alternate nucleotide. No
other base may differ. The substitutions must be at distinct, verified insert
coordinates separated by at least three bases. Strand, allele identities and
actual inserted sequences must be resolved from the source; labels alone are
insufficient. Indels and constructs with additional background differences do
not satisfy this square.

For replicate `r`, let `y00,r`, `y10,r`, `y01,r`, and `y11,r` be measured
activities on the same documented log2 RNA/DNA scale in the same cell and
matched experimental setting. The replicate contrast is

```text
E_r = y11,r − y10,r − y01,r + y00,r.
```

For literal log2 RNA/DNA ratios, this is the log2 of the ratio of the double
and reference ratios' product to the singles' product. It tests multiplicative
additivity on that scale; an additive-scale contrast answers another question.
An additive offset shared by all four activities cancels exactly. Unequal
construct-specific offsets, separate per-construct standardizations, different
cell types, and incompatible assay batches need not cancel. A source's fitted
model or transformation can also change the meaning of the contrast; its
algebra and normalization must be established before outcomes are interpreted.
Swapping either reference/alternate orientation changes the sign of `E`; a
consistent orientation must be frozen for signed results.

For the current fixed-length canonical 3mer-frequency-plus-GC linear model, no
3mer spans both edited sites at this distance, so the ideal model contrast is
exactly zero regardless of coefficients. That is a property of its representation,
not a statement about a reporter, endogenous regulation, or the genome. A
comparison of its squared interaction error with a zero-prediction baseline is
uninformative: their predictions are identical. An average *signed* contrast can
also be zero when substantial positive and negative interactions cancel.

## A possible reproducibility estimand

If a source supplies at least two genuinely independent, matched replicates per
square and cell, write `E_icr = θ_ic + ε_icr` for square `i`, cell `c`, and
replicate `r`. Here `θ_ic` is that square's mean reporter-scale interaction in
the specified experimental population. If errors have conditional mean zero and
cross-replicate covariance zero, then

```text
T_ic = mean over r < s of (E_icr × E_ics)
E[T_ic | θ_ic] = θ_ic².
```

The corresponding candidate population quantity is a prespecified weighted
mean of `θ_ic²`, for example equal weight per dependency component, equal weight
per square within a component, and equal weight over a specified set of cells.
Report its units as squared log2 activity units, and optionally its square root
as a descriptive scale. In a finite sample the estimated mean of `T_ic` can be
negative; retain it rather than truncating it to zero. A positive lower bound for
this mean, from a method justified for the actual sampling design, would support
reproducible nonzero interaction in the measured reporter system under the stated
assumptions. It does not identify a mechanism or establish an endogenous effect.

This cross-product is a candidate, not a universal correction. Independent
sequencing lanes or barcodes do not by themselves establish independent
biological replicates. A shared input-DNA estimate, construct-specific barcode
bias, shared batch effect, or shared normalization error can give nonzero
cross-replicate covariance, making `E[T_ic] = θ_ic² + Cov(ε_icr, ε_ics)` rather
than `θ_ic²`. The same problem arises from systematic construct-specific bias.
If the source cannot establish the needed independence or provide a valid
covariance model, describe observed contrasts and stop confirmatory inference on
this estimand. If true interactions vary between replicates, the target is the
square of their mean only under the corresponding zero-mean, uncorrelated
deviation assumptions, not the mean of their squared interactions.

Inference must account for squares sharing reference or single constructs,
loci, or related sequences. Define dependency components from construct
identity, genomic origin and a source-appropriate sequence-similarity rule
before looking at outcomes. Freeze the sampling unit, weights, number of
independent components and a component-level interval procedure; a component
bootstrap cannot create independent experiments or repair shared technical
errors. Sparse clusters, dependence between clusters or too few independent
replicates require a prespecified stop or a justified alternative, not a nominal
interval. Present all eligible contrasts and their uncertainties, support and
missingness along with the aggregate; do not select examples by significance.

## Before observing new outcomes

The source-specific protocol needs to bind the downloaded source and license,
design and inserted sequences, cell contexts, activity and normalization
semantics, expected replicates and covariance assumptions, square construction,
separation and orientation, prior exposure, inclusion and exclusion rules,
missing-value and low-support handling, weighting, uncertainty method, numerical
implementation, minimum support and explicit stop conditions. Pseudocounts,
source QC and published transformations must be described and fixed before
interpreting values; a missing result is not evidence of zero activity. Grouping
is done before outcome-based filtering. Independent code should reconstruct
the eligible square identities and the resulting contrasts, and preserve failed
gates and discrepancies.

Prior exposure to an assay is relevant to claims of validation. Previously used
HBB/Kircher, Agarwal and Siraj main-panel outcomes are explored evidence, even
though the exact model identity itself is outcome-free. Independent analysis of
the same source is not independent biological replication. Neither a positive
reporter contrast nor a precisely estimated near-zero contrast establishes
endogenous causal function, clinical significance, or a novel mechanism.

## Synthetic checks

Run `python3 -m unittest discover -s research/epistasis-design -p 'test_*.py'`.
The controls check offset cancellation, signed cancellation, a zero expected
cross-product with independent noise, its bias under shared error, and an
observed negative cross-product. They also check separated edits in a small
exhaustive feature enumeration and exhibit a distance-two interaction. These
synthetic checks do not validate a source, infer replicates, or replace the Lean
proof and its documented rational-to-float limitation.
