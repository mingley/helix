# Context specificity v2: frozen published-ESS evaluation

The machine authority is
[context-specificity-v2.json](context-specificity-v2.json), SHA-256
`8f0576845fecb5cfcf9bfa354fab446d849244448a3489905a6ce33bdeb21cec`.
This protocol is separate from
[blocked context prediction v1](../research/context-v1-source-gate.md).
V1 required raw measured activity; the permitted files instead supply design
information and derived specificity scores. V1 remains unchanged and has no
model fit or activity-prediction result.

## Exposure and target

Both workbooks were acquired under v1. Their headers, IDs, design DNA and
coordinates were inspected before this specification. The protocol author did
not inspect numeric ESS values or fit a model. This is an outcome-unseen
analysis specification following a source-availability failure; it is not an
untouched payload or a prospective experiment.

The source is Agarwal et al., DOI `10.1038/s41586-024-08430-9`, with the
preprint-associated [Zenodo 10558183 release](https://zenodo.org/records/10558183).
It declares CC BY 4.0. The source panel was deliberately selected using prior
activity bins and variable-promoter criteria. The source paper already studies
and predicts ESS. This study tests a small fixed model under stronger
loci/homology grouping; it does not discover that regulatory context matters.

The target is the **published Element Specificity Score**, unchanged, in HepG2,
K562 and WTC11. Under the paper's “Calculation of element specificity scores”
methods, per-cell log2 RNA/DNA activities were z-score-transformed and each
element's across-cell mean was subtracted:

`ESS_ic = z_ic - (z_iHepG2 + z_iK562 + z_iWTC11)/3`.

The subtraction uses all three cells, not the other two. The cited methods do
not establish an independently reconstructable source SD convention or exact
normalization population size. Use the published values; do not reconstruct
raw activities, source z-scores, source means or source SDs.

**The source's panel-wide target normalization includes information from future
test elements.** Our grouped withholding controls model fitting, not that prior
target construction. No prospective or leakage-free outcome-preprocessing claim
is made. ESS is dimensionless specificity, not raw reporter activity or a
log2-effect coefficient. There is no additional target normalization.

## Exact assets and schema

| Workbook | SHA-256 | Published MD5 |
| --- | --- | --- |
| Table S7 - joint lib design.xlsx | `f9ec4c6c2d230b0ef507a39e5d5c73b3b36f90ccc5758c33b442c7561e20e7e5` | `5d9b71452e9648e8473edd74a4ae26fc` |
| Table S9 - folds and performance.xlsx | `0f9ff22d4b46fab36d419399e6ab743dc08a0813f812dc7bd96f90eb7e0f2b8a` | `ef967c8388d84da9ec81d22e11c868bd` |

Verify these bytes before parsing scores. Only these two workbooks supply
design/labels. CNN archive predictions, `prediction.N` and recomputed surrogate
labels are ineligible. Authoritative hg38 reference acquisition is permitted
for design validation with a pinned source/checksum receipt.

S7 worksheet `joint library` uses header row 2. Join its exact `name` to S9's
exact `Sequence ID` in worksheet `Fold assignment, joint library`, header row 1.
Use exactly these S9 score columns:

- `HepG2, Element Specificity Score`
- `K562, Element Specificity Score`
- `WTC11, Element Specificity Score`

The S9 `Fold` column is ignored for our splits. No alias search or suffix/fuzzy
ID transform is permitted. A pre-outcome ID-only audit found 60000 unique S7
IDs and 55338 unique, nonempty S9 IDs, all exact unique S7 matches. No `_F`
fallback was needed. S9's maximum worksheet dimension implies 59107 possible
data row numbers but is not a record count: its XML has 55338 nonempty ID cells
and one explicit blank-ID row after the header. Retain gaps and blank rows in
the audit. The ID-audit receipt SHA-256 is
`1035d86ab6eaa97d01dd829f4298fcefb98b4a65eb2bc0a92f4b536bc6402ac2`.

The explicit S7 sequence header is
`230nt sequence (15nt 5' adaptor - 200nt element - 15nt 3' adaptor)`.
Require full length 230 and extract the core using zero-based slice `[15:215]`.
Record both 15-nt adapters and full-construct/core hashes; features and reference
validation use only the 200-nt core.

Use `chr.hg38`, `start.hg38`, `stop.hg38` and `str.hg38` as independently
verified zero-based half-open coordinates. Require nonmissing coordinates,
span exactly 200 and an exact core match to the pinned hg38 reference in the
declared strand, reverse-complementing reference for a minus-strand insert.
All verified declared hg38 contigs, including alternate contigs, are eligible.
Record any one-to-one contig-name alias. No silent coordinate offset, strand
correction, liftover-gap rescue or DNA replacement is allowed.

Design-only inspection found 59319 span-200 intervals, 11 non-200 intervals and
670 missing-coordinate designs among 60000 records. Non-200 and missing
intervals fail the fixed reference gate; no favorable rescue is attempted.
Retain all source-only, invalid/reference-mismatched, missing/nonfinite and
incomplete-context rows with reasons. No category, ESS magnitude/rank, source
fold or zero-sum-residual filter is allowed. Source replicates/barcodes are not
additional elements. Exact/RC duplicate source IDs remain but share a group.

## Fixed grouping, folds and models

The grouping, fits and primary loss weighting follow v1's fixed design, with
the now explicit 200-nt DNA core and ESS target:

- Union documented parents, exact/RC duplicates, positive-length genomic
  overlaps, all touched 1-Mb bins anchored at zero and equal-length pairs with
  at most `floor(L/10)` Hamming mismatches in either orientation. S7 has no
  separate parent field; do not infer a family from a broad cell/category prefix.
- Use exhaustive pigeonhole candidate generation, exact verification and
  transitive components. Abort above 100 million unordered candidate pairs;
  never drop large buckets or continue an incomplete grouping.
- Require at least 100 components and largest component at most 40%. Assign
  five folds greedily by decreasing component size then SHA-256 of sorted IDs,
  choosing the least-populated fold and breaking ties by fold index.
- Exclude entire training components within 10 kb of any test interval. Each
  guarded training fold must retain at least 50% of all elements and all three
  contexts. Keep every test element. Any feasibility failure aborts the primary
  study rather than weakening the rule after seeing outcomes.

Features are 32 lexicographically ordered RC-canonical trinucleotide
frequencies followed by GC fraction, using only core DNA. Count overlapping
3-mers and divide once by 198; GC is `(G+C)/200`. One training-only population-SD
transform is shared across contexts; SD at most `1e-12` becomes scale one.
Source categories, activities, folds, predictions, counts and coordinates are
not predictors.

Let `m` be the grand training mean, `m_c` each context's training mean and `z`
standardized DNA. Use five fixed models:

| Model | Prediction |
| --- | --- |
| G: global mean | `m` |
| C: context means | `m_c` |
| S: shared DNA | `m + z beta` |
| A: additive DNA/context | `m_c + z beta`, sharing S's stored slopes |
| I: independent-context DNA | `m_c + z beta_c` |

S fits the `n` element-wise context-mean targets with SSE plus
`10 ||beta||²`, equivalent to a pooled penalty of `3*10`. Each I fit uses its
`n` context-specific targets and penalty ten. Intercepts are unpenalized; no
tuning, feature selection or response rescaling is allowed. A is computed
directly as context mean plus the shared dot product.

Source centering implies each element's ESS mean should be approximately zero.
Consequently G/S should be near zero and C/A nearly identical. **Fit the actual
values anyway:** do not project ESS to zero sum, snap fitted slopes/intercepts
or discard residual-violating rows. Report signed/absolute/RMS/maximum row-sum
residuals and counts above each fixed threshold `1e-12`, `1e-8`, `1e-6`, along
with fit magnitudes and rounding provenance. Investigate semantic contradictions
openly. Also report an additional fixed prediction of zero as an audit baseline;
it is not a selected sixth fitted model.

## Primary comparison and uncertainty

Primary MSE averages the three context errors within each assayed element and
then averages elements equally. MSE units are ESS squared; MAE uses ESS units.
Components determine dependence/splits, not primary equal-group reweighting.
Report every model/context loss and equal-component MSE as a fixed sensitivity.

The **single primary contrast** is `MSE(A)-MSE(I)`, positive favoring context
slopes beyond additive offsets. The retained `MSE(S)-MSE(A)` and
`MSE(C)-MSE(A)` contrasts are fixed diagnostics, not co-primary tests: source
normalization already constrains their targets. Retain all results, favorable
or unfavorable, and the fixed-zero losses.

Use 2000 common component-bootstrap draws, SplitMix64 hexadecimal seed
`0x2026102702` (decimal `138077546242`). Sample `K` component IDs with replacement
in ascending component-ID order. Pool member loss sums and counts, with draw
multiplicity, so each draw remains element-weighted. Reuse the same draws for
all three contrasts. Type-7 2.5th/97.5th percentiles are 95% conditional ranges
for fixed out-of-fold predictions. They omit refitting, target-normalization,
source-selection and study-transfer uncertainty. No calibrated significance or
p-values are claimed.

Report all three observed/predicted cell-pair differences with A/I MSE and MAE.
Ridge linearity on the balanced design implies mean I slopes equal S slopes,
and mean-context I/A predictions agree. Thus, at identical element weights,

`MSE(A)-MSE(I) = mean_pair[contrastMSE(A)-contrastMSE(I)] / 3`.

This is a structural consistency diagnostic, not independent corroboration.

## Verification and limits

Independently reproduce input/schema binding, eligible IDs, integer feature
counts/lengths, groups, folds, guards and bootstrap draws exactly. Prediction
absolute tolerance is `1e-7`, loss tolerance `1e-8`. Preserve numerical failures;
do not widen tolerances or require exact-rank agreement. Perturbing test-only
outcomes in synthetic tests must leave every fitted quantity unchanged.

Publish the full source/exclusion/grouping/attempt ledger, reference receipts,
all fixed models and diagnostics, primary contrast and independent-comparison
status. Positive predictive performance would concern this small model,
selected panel and fixed, source-normalized specificity target in known cells.
It does not recover raw activity, identify a causal mechanism, establish new
biology or support unseen-cell or independent-study generalization. Any later
change requires a separately versioned protocol and exposure disclosure.
