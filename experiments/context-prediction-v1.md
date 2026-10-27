# Context prediction v1: frozen acquisition and analysis protocol

The machine authority is [context-prediction-v1.json](context-prediction-v1.json),
SHA-256 `a089a6668ad73e943797df405b3ba510b1377190b416da1bd8df12bc368b77c2`.
The analysis was specified before acquiring the listed source payloads. This
is a retrospective public-data study of a deliberately selected panel, not a
prospective experiment, raw-count preprocessing pipeline or independent
validation of the source study. Its original [design draft](../research/context-evaluation-design.md)
records the reasoning; this frozen JSON resolves the draft's implementation
choices.

## Scientific question and source gate

The question is whether known assay context improves activity prediction for
held-out sequence families beyond context means and shared DNA composition, and
whether context-specific DNA slopes improve beyond additive context offsets.
It concerns new elements in three observed contexts, HepG2, K562 and WTC11. It
does not test an unseen context or an independent source study.

The source is Agarwal et al., DOI `10.1038/s41586-024-08430-9`, with the
preprint-associated Zenodo release [10558183](https://zenodo.org/records/10558183).
The repository declares CC BY 4.0. Preserve the actual license/attribution
snapshot and stop if score-data reuse is contradicted. Prior activity-bin and
variable-promoter selection restrict inference to the selected panel. The
source already studies context-dependent prediction; that general concept is
not a novel result of this analysis.

The only permitted source payloads are:

| Asset | Bytes | Published MD5 |
| --- | ---: | --- |
| Table S7 - joint lib design.xlsx | 8301566 | `5d9b71452e9648e8473edd74a4ae26fc` |
| Table S9 - folds and performance.xlsx | 3564742 | `ef967c8388d84da9ec81d22e11c868bd` |
| 30 explicitly listed joint-result TSV members | Member sizes/CRC in JSON | Whole ZIP: `edf24033ce6adf8e3cb27e465730c329` |

Inspect the bounded workbooks first, S7 then S9. The 30 TSVs are candidate
provenance sources, not assumed observed labels: source model output may
contain only IDs and predictions. The 3967659598-byte ZIP has not been
downloaded whole. Partial range acquisition must verify local/central metadata,
decompression, member sizes and CRC32 and retain HTTP/range receipts and member
SHA-256 hashes. It cannot claim to verify the whole-archive MD5. The frozen
central-directory index SHA-256 is
`42e5032c42a6abfc79a1574dae2f10cf17b4bb094a183a66412e948e8a69522c`.

**Fail closed:** the label must be observed activity, specifically the arithmetic
mean of three published per-replicate median-normalized log2 RNA/DNA activities,
or a documented published aggregate exactly equal to it. Source CNN predictions,
`prediction.N`, ensembles, ESS, full-panel z-scores and design activity bins are
never empirical labels. If the listed files cannot unambiguously supply this
endpoint, actual assayed DNA and its coordinates, publish a blocked-acquisition
report and stop. Do not silently substitute a score, cohort or source.

Headers and authoritative methods may bind the fixed semantic contract to exact
columns after acquisition. Save a source-hashed schema-binding artifact before
fitting or outcome-pattern inspection. A column name alone does not prove its
meaning. Binding may instantiate fields and documented missing-value tokens;
it cannot change the endpoint, transform, cohort source or analysis. Conflicting
legitimate measured copies require resolution or a stop, rather than choosing
the more favorable values.

Table S7 is a preprint-associated table, not automatically final-paper
Supplementary Table 10. Establish version identity explicitly. Published median
normalization may use the whole source experiment; our training-only model
processing does not make that source normalization prospective. No additional
all-panel z-score, response centering or raw-count reconstruction is allowed.

## Cohort, groups and folds

Use all source joint-panel element IDs with the same actual assayed insert,
valid uppercase A/C/G/T DNA, unambiguous verified same-assembly coordinates and
finite required observed activity in all three contexts. A common parent label
with different inserts is insufficient. Keep missing/incomplete/invalid records
and their reasons in a separate audit; do not impute or filter on activity,
p-value, replicate agreement or predictability. Replicates and barcodes are not
independent elements. Distinct exact/RC duplicate source IDs remain in the
panel but share a group, and their multiplicity is reported.

Build connected components by unioning documented parents/families, exact or
reverse-complement DNA identities, positive-length genomic overlaps, all touched
1-Mb bins anchored at zero, and equal-length pairs with at most `floor(L/10)`
Hamming mismatches in either orientation. The pigeonhole index in the JSON is
exhaustive for that last rule. Abort above 100 million distinct unordered
candidate pairs; do not skip low-complexity buckets or continue a partial
grouping. Unequal-length or indel homology is not guaranteed to be captured.

Require at least 100 components and no component exceeding 40% of the primary
panel. Hash each component's newline-joined sorted IDs with SHA-256. Sort by
decreasing element count then hash, and greedily place it in the least-populated
of five folds, with fold index breaking ties. All contexts stay together.

For each held-out fold, remove a whole training component if any member is
within 10 kb of any test interval on the same chromosome. Half-open interval
distance is zero for overlap or adjacency. Every guarded training fold must
retain at least 50% of all primary elements and all contexts. Failure aborts
the primary study; thresholds are not relaxed after outcomes. Save complete
group/fold/guard ledgers. Source testfold names are not this study's folds.

## Five fixed models

Features are 32 lexicographically ordered reverse-complement-canonical
trinucleotide frequencies followed by GC fraction. Count overlapping 3-mers in
integers and divide by `L-2`; GC is `(G+C)/L`. Use only actual insert DNA.
Lengths, source predictions, counts, tags, activity bins and outcomes are not
predictors. On each fold calculate one training-element feature mean and
population SD; scales at most `1e-12` become one. Apply that same transform in
every context. No tuning or feature selection is allowed.

Let `m` be the grand training mean, `m_c` the training context mean, `z_s` the
standardized DNA features and `a_s` the mean observed training activity across
contexts. All intercepts are unpenalized.

| Model | Prediction |
| --- | --- |
| Global mean G | `m` |
| Context mean C | `m_c` |
| Shared DNA S | `m + z_s beta` |
| Additive context A | `m_c + z_s beta`, sharing S's stored slopes |
| Independent context I | `m_c + z_s beta_c` |

Fit shared slopes on the `n` sequence-mean targets using sum-squared error plus
`10 ||beta||²`; a pooled `n*C` representation would require penalty `C*10`.
Fit each context's `n` targets with penalty ten. Both use
`(Z'Z+10I) beta = Z' centered_target`. Responses remain in native log2 units.
Compute A directly as context mean plus the shared dot product.

The balanced panel implies S and A have identical slopes and the mean I slope
equals S's slope. Mean-context I and A predictions must therefore agree within
the fixed numerical tolerance. These are predeclared artifact checks, not
outcome-based model changes. I has more parameters; a positive held-out result
would support this composition model's context dependence, not identify a
causal molecular mechanism.

## Scoring and conditional uncertainty

Primary MSE averages the three context errors within each element, then averages
elements equally. Components determine dependence and splits; they do not make
the primary an equal-component mean. Report MAE, every context's losses, all
five models and equal-component MSE as a fixed sensitivity.

Prespecified improvements are baseline MSE minus candidate MSE:

1. S minus A: context beyond shared DNA.
2. C minus A: DNA beyond context means.
3. A minus I: sequence-by-context slopes beyond additive offsets.

Keep all three contrasts, including unfavorable results. Context offsets alone
are insufficient evidence of a context-dependent DNA relationship.

Use 2000 common component-bootstrap draws with SplitMix64 seed
**`0x2026102702` = `138077546242`**. This is a hexadecimal seed, not the earlier
regulatory study's decimal seed. Order components by ID, sample `K` component
indices with replacement for each draw and use identical draws for every
contrast. Pool full member loss sums and counts, with multiplicity, to retain
equal-element weighting. The first three unsigned outputs are
`311730414543f7d5`, `b0ee1155b9089752`, `fea23a01dc5344f1`.

Report Type-7 2.5th/97.5th percentiles as 95% **conditional resampling ranges**
for the fixed out-of-fold predictions. They omit refitting uncertainty,
residual dependence, source selection and study-to-study variation. They are
not calibrated significance tests or genome-population confidence guarantees.
No null p-values or significance declaration is part of this version.

For all three context pairs report MSE/MAE of predicted versus observed activity
differences for A and I. With three complete contexts and matching weights,

`MSE(A)-MSE(I) = mean_pair[contrastMSE(A)-contrastMSE(I)] / 3`.

This identity isolates the context-varying component because mean-context
predictions agree. All pairs are retained; it is a consistency check and
description of the interaction, not another independent experiment.

## Reproduction and failure reporting

Independently reconstruct inputs, schema binding, eligibility, integer feature
counts/lengths, groups, folds, guards and bootstrap draws exactly. Independent
predictions must agree within absolute `1e-7`, loss metrics within `1e-8`.
No exact-rank endpoint is required. Synthetic tests must also show that changing
test-only outcomes cannot affect fitted quantities and that the balanced-model
identities hold.

Retain all acquisition failures, aborted grouping, missingness, failed numerical
checks, attempts and negative contrasts. A failed semantic source gate is a
valid documented outcome of this protocol. Any substituted file, endpoint or
analysis requires a separately versioned protocol and exposure disclosure; the
announced JSON bytes remain unchanged.
