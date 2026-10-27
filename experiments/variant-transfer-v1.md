# Frozen variant-effect transfer protocol

The [machine-readable protocol](variant-transfer-v1.json) tests whether the fixed
Agarwal activity model predicts signed single-nucleotide effects in a separate
Siraj reporter experiment. This is a repository freeze, not external
preregistration. Its commit must precede the first S4 acquisition request.

Protocol SHA-256: `811930c2d3a5ee763014ea611acda4d2c84404b7d5ae399c7843c62d0a828149`. Frozen at
`2026-10-27T22:21:30.598229+00:00`.

## Source and prior exposure

Siraj et al., *Functional dissection of complex trait variants at
single-nucleotide resolution*, Nature (2026), DOI
[10.1038/s41586-026-10121-6](https://doi.org/10.1038/s41586-026-10121-6).
The selected sources are publisher-linked S2 variant oligo sequences and S4
main-panel MPRA results. This protocol excludes the saturation-mutagenesis
posterior effects, multivariate-shrunken estimates, source model predictions,
diplotypes and across-cell meta-analysis.

Earlier HBB/Kircher outcomes and Agarwal specificity/activity values and fitted
models are known. Siraj methods and aggregate metadata have been read. DNA-only
S2 was acquired before this freeze: **40,795,798 bytes**, SHA-256
`4b6e20d1789021ee8f11f57af366a65680d25496de9351a9e781e9424716d824`.
Its sheet `S2.Variant Oligo Seqs` has title row 1 and header row 2:
`Variant`, `Allele 1 Oligo`, `Allele 2 Oligo`. The source reports 365,762 data
rows; that is not the final eligible cohort. Before the final semantic freeze, DNA-only inventory found 336,475 exact
200 nt primary-chromosome SNV pairs and 16 prior-panel exact/RC allele-match
records. These counts do not define the admitted cohort. Source-build/reference
validation and full prior-source overlap checks remain pending. **No numeric S4 outcome payload was acquired or
parsed at this freeze.**

The exact official S4 URL is bound in the JSON. Its literal sheet/header mapping
is not yet known. Expected producer names such as `Log2Skew` or `log2Skew` are
schema expectations, not an authenticated spreadsheet binding. Source-selected
loci and overlapping investigators limit independence; this is a separate
experiment, not a random genome sample or independent proof of a mechanism.

## Required acquisition and second freeze

1. Commit this semantic protocol before the first S4 GET. Preserve normal proxy,
   CA trust and TLS verification, with finite connection/read/total timeouts.
   Download at most **400,000,000 compressed bytes**; allow at most **10,000 ZIP
   members**, **2,000,000,000 total XML bytes** and **2,000,000,000 total
   uncompressed bytes**. Reject HTTP/HTML errors, wrong magic, incomplete bodies,
   encryption, duplicate names, unsafe paths and symbolic links. Check every
   member CRC and actual size within those limits. CRC streaming does not
   authorize interpretation of numeric response cells.
2. Inspect archive metadata, textual headers and source legends/license only.
   Retain the exact file hash and acquisition receipt. Do not read activity or
   allelic-effect values to choose a column.
3. Bind exact IDs, K562/HepG2 labels, allele orientation, units, header positions,
   DNA/adaptor fields, source genome build, QC, duplicate-library choice and
   missing tokens. Freeze the implementation, metadata-eligible ID list,
   exclusions, complete explored HBB/Kircher inventory, groups, source/header
   fingerprints and all five old fit hashes in a **second committed receipt**.
   Bound coordinates, source classes and support metadata may be read for these
   checks; response cells remain unopened.
4. Only then parse the selected response columns and execute the fixed comparison.
   Record documented missingness and every failure. A source-semantic ambiguity
   blocks this protocol.

The article states CC BY 4.0; confirm that the selected workbooks have no contrary
credit/terms before redistribution. The separate Zenodo deposit states
CC BY-NC 4.0, so its terms are not assumed to match the article. Keep original
large source bytes outside Git when necessary, with exact acquisition/hash
receipts and independently readable derived bindings.

Textual bindings may implement the existing semantic choices. Any semantic
change to source, endpoint, orientation, contexts, length, eligibility, QC,
grouping, guards, model, metric, uncertainty, tolerance or acquisition bounds
requires a new explicitly exposed version **before any new numeric outcome
parsing**. This frozen file is not edited to accommodate a convenient source.

## Endpoint and admitted population

Use the source's cell-specific DESeq2 interaction estimate of **Alt minus Ref
RNA/DNA activity**, in published summit-shift-normalized delta log2 fold-change
units. Preserve literal source floats. There is no target centering, scaling,
clipping, effect-oriented sign recoding or calibration to the training study.
The primary includes **both K562 and HepG2**; neither cell can be selected after
observing its results.

Eligible pairs must be human primary-chromosome SNVs with exactly one nucleotide
change between two exact **200 nt**, uppercase ACGT genomic cores. Remove only
explicitly documented 15 nt adaptors if present. Exclude indels, ambiguous bases,
padded/inferred cores, synthetic/shuffled controls and previously reused
technical controls. Include genuine fine-mapped, eQTL/trait and human
location/annotation/null-control SNVs without a PIP, CRE, activity, emVar, FDR or
effect-magnitude filter.

Respect published DNA-count/RNA-count QC and source coverage-based duplicate
library selection. The source methods exclude either allele with DNA counts
below 30 or RNA counts equal to zero. An already prefiltered, deduplicated S4
must be used unchanged and described as source preprocessing. Do not reconstruct
unknown counts or choose a library by activity. Conflicting or unexplained
duplicates block evaluation.

The primary complete-case cohort requires ordinary finite numeric effects in
both declared cells. Only source-documented missing tokens may exclude an
otherwise eligible variant; retain every ID and missingness count by cell and
component. Formulas, errors, booleans, numeric strings, unexplained missingness
and contradictory IDs are gate failures. At least **100 complete pairs and
10 components** are required for the primary interval claim; otherwise retain
descriptive results and report the support failure.

## Whole-component isolation

All observations are external test data; no new-study training split is needed.
Before responses are parsed, join variants transitively by exact/RC DNA identity,
source parent/locus, overlapping intervals, **1 Mb bins** and **at least 90% full
200 nt ungapped identity in either strand**, considering both alleles. Key bins
by reference build/chromosome and zero-based interval-start integer division by
1,000,000. Shared loci/homology may join bins; preserve those components intact.

Exclude an entire new component if any member overlaps or lies within the
inclusive **10,000 bp** interval guard of any old Agarwal training locus, or
matches an old training core by the declared identity/homology rule. Consider
the union of all five old training sets. Apply equivalent exclusion to a complete
explicit inventory of all explored HBB/Kircher loci and sequences; an incomplete
inventory blocks the study. Bind genome build and coordinate convention and
verify converted reference cores exactly. Do not compare hg19 coordinates
directly with hg38. Missing/ambiguous mappings are recorded before responses.

Freeze the exhaustive/indexed candidate implementation and budget in the second
receipt. Do not omit difficult homology candidates to finish a run.

## Fixed predictors and primary comparison

Use the five published independent-context activity fits from
[the clean activity core](../research/results/context-activity-v3/core/artifacts.json),
artifact-manifest SHA-256
`2ea476bc5413dc13486dec9f3d7e9162a50a0e47f1d83ce29efabc5f5cc5af9d`,
from code commit `035ca5731fe594622cd2c3c27cd3bb51ec65c924`.
Their source protocol SHA-256 is
`64e14037fb1f42ddc1ecaca1d23bb8096fabac1f1e892d105eaa20f525fbbb20`.
No new outcome model or GC-only model is fitted.

For each cell, average the five fixed model differences:

```text
predicted_delta = mean_fold sum_feature (
    standardized_coefficient / training_feature_scale
    * (feature_alt - feature_ref)
)
```

The existing 32 RC-canonical overlapping 3mer frequencies and GC fraction retain
their exact order; denominators are 198 and 200. The old intercepts and feature
means cancel. There is no new coefficient fitting, penalty selection or cross-study
response calibration. Shared A/S slopes supply one declared secondary predictor.

The **single primary comparison** is fixed matched-cell I prediction against
**zero effect**. Average squared errors equally across K562/HepG2 within a
variant, equally across variants within a component, then equally across
components. Positive macro MSE improvement means the candidate beats zero.
Evidence of useful transfer requires a positive observed improvement and a
positive lower endpoint of the declared conditional interval. A failed primary
cannot be replaced by a favorable secondary metric or subgroup.

## Uncertainty, controls and limits

Use **2,000 paired whole-component bootstrap draws**, sampling the lexicographically
sorted component list with replacement using verified SplitMix64 seed
`0x2026102801`. Carry every member/context of each sampled component. Report
Type7 2.5% and 97.5% quantiles and retain all indices/statistics. These are
fixed-predictor conditional selected-panel ranges, not calibrated p-values or
biological replication.

Secondary results include MAE, per-cell comparisons, equal-variant sensitivity,
shared-slope prediction, Pearson signed correlation, exact sign agreement with
zeros retained, and metadata-defined GC-preserving/changing SNVs. No significance
or near-zero filter is allowed. Synthetic controls check identity, allele-swap
antisymmetry, RC invariance, exact local trimer changes and zero predicted
interaction between substitutions separated by at least three positions.

An independent implementation must read the original XLSX, rederive eligibility,
DNA features, fixed predictions, component losses and resampling. IDs, integer
counts, source float bits and resampling indices must agree exactly; prediction
and metric maximum absolute tolerance is **1e-10**, without a relative term.
Keep all failed checks and earlier attempts; never widen tolerances to pass.

This evaluates transfer of a small composition model on selected reporter
constructs. It does not establish endogenous regulatory causality, clinical risk,
absolute expression, a universal mechanism or global scientific novelty.
