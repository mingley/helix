# A replicate-aware test of equal-feature construct predictions

This is a **source-agnostic methods note**, not a frozen analysis protocol, a
claim of a biological effect, or a release of reporter outcomes. It describes a
candidate estimand for the feature-identical constructs found in the
[S18 design-only audit](../epistasis-feature-collisions-v1/README.md). A
source-specific protocol must establish the actual measurements, replicate
structure, normalization, selection, and uncertainty procedure before applying
these equations. The Siraj study and other outcomes from it have been explored
in this project; this design-only calculation does not make that study wholly
unexposed or establish that any matching replicate values are available.

## What a collision forces

Fix an assay context and a complete model input. It must include any cell,
library, orientation, batch or other information that a proposed predictor is
allowed to use. Within that context, a feature class contains distinct DNA
constructs whose complete inputs are identical. Every deterministic predictor
using only those inputs returns one common prediction `t` for that class,
whether the predictor is linear or nonlinear. A model given additional sequence
or design information may distinguish its members. Results in one cell or
assay do not establish equality or inequality in another.

For construct `j`, let `theta_j` be its mean activity in a **specified**
population of genuinely repeated experimental units, on one verified common
scale. In this note it is an assay-specific mean, not endogenous activity or a
clinical effect. Let `w_j >= 0` be fixed construct weights within the class,
with `sum_j w_j = 1`, and let `theta_bar = sum_j w_j theta_j`. Then

```text
sum_j w_j (t - theta_j)^2
  = (t - theta_bar)^2 + V
  >= V,

V = sum_j w_j (theta_j - theta_bar)^2
  = sum_{j<k} w_j w_k (theta_j - theta_k)^2.
```

`V` is the minimum weighted squared error over an unrestricted common scalar
prediction for that class; a restricted model may have greater error. For two
constructs with equal weights it is `(theta_1 - theta_2)^2 / 4`. This is a
mathematical lower bound on the stated *latent-mean* risk for the stated finite
class and weights. It is not, by itself, a bound on endogenous function or on
generalization to other sequences. If observed labels are used directly in the
same identity, the resulting bound includes their measurement error and is not
automatically the latent-mean bound.

The pair form is useful computationally, but weights matter. For four equally
weighted members, `V = sum_{j<k}(theta_j-theta_k)^2 / 16`; simply adding the
six separate two-construct bounds would quadruple this class variance. Classes
and their members must be defined without duplicating the same observation
unless a specified weighted multiset is deliberately the target. A weighted
average of separate pairwise bounds is a bound for its induced endpoint weights;
it is not automatically a bound for a differently weighted dataset.

## Cross-replicate estimator

For a pair in a class, suppose `R >= 2` matched replicate units each measure
both constructs in the same cell, setting and compatible activity units. Write

```text
D_r = Y_{j,r} - Y_{k,r} = Delta + epsilon_r,
Delta = theta_j - theta_k.

U_{jk} = mean_{r<s}(D_r D_s)
       = ((sum_r D_r)^2 - sum_r D_r^2) / (R (R - 1)).
```

If every `epsilon_r` has conditional mean zero and each distinct pair of
replicates has conditional covariance zero, then
`E[U_{jk} | theta] = Delta^2`. With fixed class weights,
`V_hat = sum_{j<k} w_j w_k U_{jk}` is then unbiased for `V`. A common
additive offset applied to *both constructs in the same replicate* cancels in
`D_r`. Construct errors may be correlated **within** a replicate, and variances
may differ by construct or replicate; neither fact alone breaks this expectation.
The cross-replicate mean and covariance assumptions still need to hold for the
pair difference. Finite moments and an appropriate sampling model are needed
for uncertainty calculations.

This estimator can be negative in a finite sample, including for `R = 2` when
the two observed differences have opposite signs. Retain negative estimates;
do not replace them with zero. An unbiased estimate of `V` is not itself a
deterministic lower bound on the unknown realized risk in a particular sample.
The estimate is in squared units of the specified activity scale. Taking a
square root changes its statistical properties and is not a generally unbiased
estimate of a root-mean-square difference.

If replicate errors are centered but cross-replicate covariance is not zero,
the precise bias is

```text
E[U_{jk} | theta] = Delta^2 + mean_{r<s} Cov(epsilon_r, epsilon_s | theta).

Cov(epsilon_r, epsilon_s)
 = Cov(e_{j,r}, e_{j,s}) + Cov(e_{k,r}, e_{k,s})
   - Cov(e_{j,r}, e_{k,s}) - Cov(e_{k,r}, e_{j,s}).
```

Here `e_{j,r}` and `e_{k,r}` are the construct-level errors after the stated
means. This covariance can be positive or negative. Nonzero conditional mean
errors add further bias: if `m_r = E[epsilon_r | theta]`, a cross-product's
expectation contains `Delta (m_r + m_s) + m_r m_s` as well as covariance.
Sequence-dependent systematic bias cannot be assumed to disappear merely
because there are multiple observations.

As an example, suppose both constructs use one shared input-DNA preparation and
`Y_{j,r} = theta_j + U_j + e_{j,r}`. The common construct-specific difference
`U = U_j - U_k` remains in every `D_r`. If `U` is centered conditional on
`theta` and independent of centered, cross-replicate-uncorrelated `e`, the
expectation over new input preparations is `Delta^2 + Var(U | theta)`, not
`Delta^2`; in a particular retained preparation the conditional target is
`(Delta + U)^2`. If the scientific target is explicitly conditional on that
fixed preparation, its conditional mean may be a legitimate descriptive target,
but it does not establish reproducibility over new DNA preparations. A shared
error that affects both constructs equally cancels; a shared denominator,
barcode composition or normalization does not in general guarantee that
equal effect. Sequencing lanes or barcodes from one library are not, by
themselves, independent biological replicates. Replicate transfections sharing a
plasmid pool may still have this shared-input covariance.

## Aggregate target and uncertainty

Define classes, contexts and eligible replicate units from design and source
metadata before inspecting selected responses. A proposed aggregate is

```text
V_target = sum_g lambda_g V_g,       V_hat_target = sum_g lambda_g V_hat_g,
```

where nonnegative `lambda_g` sum to one and all within-class and between-class
weights are fixed in advance. For example, one could give prespecified
dependency components equal total weight, then eligible classes within each
component equal weight, with explicitly chosen construct weights inside a
class. This is one possible target, not a universal weighting rule. It differs
from equal weight per square, per physical construct, per sequence, or per
genomic locus. Overlapping windows, repeated source rows, the same DNA sequence,
related loci and shared measurements can otherwise give one observation
multiple weights. Publish the complete class and weight ledger.
Singleton classes have `V_g = 0`; including or excluding them changes the
aggregate weights and target. A floor conditional on collision-enriched classes
must not be presented as the same weighted floor for the whole design population.

Replicate dependence within a pair affects the expectation above; dependence
between pairs, loci, classes and cells additionally affects uncertainty. Define
dependency components by documented construct identity, physical library,
genomic origin and a source-appropriate homology rule, with transitive links,
before outcomes. A component is not automatically an independently sampled
experimental unit. The source-specific protocol must say which random units
and populations the interval covers, how shared transfections and input-DNA
preparations enter its covariance or resampling scheme, what number of genuinely
independent units is required, and when it stops for sparse or dependent data.
For example, resampling loci conditional on a single shared DNA preparation
cannot quantify variation between independent DNA preparations. A bootstrap or
nominal interval is not justified merely because many overlapping windows exist.

A candidate confirmatory rule is to designate **one** aggregate `V_target` as
primary and require a justified, predeclared one-sided 95% lower confidence
bound strictly above zero. The interval method must have its stated coverage at
the boundary `V_target = 0` under the actual dependence and selection design;
neither the cross-product formula nor a large count of windows establishes that
coverage. A source-specific protocol must fix the method, weights, minimum
independent support and failure rules in advance. If it cannot justify an
interval, stop confirmatory inference rather than substitute a nominal test.
Other cells, subsets, metrics or thresholds must be designated exploratory or
included in a prespecified multiplicity procedure.

A positive lower bound would support nonzero between-construct differences
under the specified reporter assay and assumptions; it would not identify a
transcription factor, causal mechanism, or endogenous effect. A bound including
zero is not proof of equality. A precision or equivalence question requires its
own margin and interval rule frozen before the selected responses.

## Source gates before observing selected responses

The protocol must bind the original inserted sequences and exact complete model
inputs; cell type and assay; source files and hashes; permitted source licenses;
activity units and allele direction; per-replicate values and their provenance;
what was independently cultured, transfected and sequenced; DNA preparation,
barcodes and normalization; known prior exposure; complete class construction;
and all dependence and weighting rules. Literal source blanks, missing tokens,
zero counts, low-support constructs, partial replicates, duplicate identities
and nonnumeric values need explicit handling. A fitted or pooled fold-change
and standard error alone do not supply the matched replicate differences or
their cross-replicate covariance.

Predefine replicate matching and how a pair with fewer than two qualifying
independent units is reported. Never silently turn a missing measurement into
zero or drop a member and call its original class unchanged. Complete-case
analysis changes the target to the retained data unless additional assumptions
are justified. Missingness or quality filters based on RNA counts, measured
activity or source significance can depend on the outcome and on measurement
errors; even a seemingly technical filter can bias cross-products after
selection. Freeze the operational rule without inferring undocumented source
missingness, report attrition and support by class and component, and stop
confirmatory interpretation if the source cannot support the required mean,
covariance and selection assumptions. If available cases use different
replicate sets, their unbiasedness and weights need source-specific justification
rather than an automatic application of the formula.

The feature-collision selection is outcome-free but is already documented after
the surrounding study and its prior exposure. It does not become a new
outcome-naive laboratory experiment by writing a later protocol. A future
source-specific plan must explicitly record all exposure and distinguish
discovery, exploratory reuse and genuinely independent measurements.

## Synthetic controls

Run `python3 -m unittest discover -s research/epistasis-collision-design -p 'test_*.py'`.
Exact finite controls check the weighted class identity, the danger of summing
overlapping pair floors, expectation with independent replicate errors,
positive and negative cross-replicate covariance, cancellation of a genuinely
common offset, a shared-input counterexample, a negative realized estimate and
outcome-dependent selection bias. These checks validate the stated arithmetic,
not any assay's replicate independence, missingness or biological interpretation.
