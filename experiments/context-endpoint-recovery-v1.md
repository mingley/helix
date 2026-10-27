# Final-paper joint-panel endpoint recovery

The frozen authority is
[context-endpoint-recovery-v1.json](context-endpoint-recovery-v1.json), SHA-256
`00175e41cf9fa0cbf177ce8fc7ef143eb508aec8b11d058f7acb4458284dc868`.
The parent reviewed and committed this contract as `2996b50` before starting
the official attachment download. This human explanation was completed while
acquisition ran; the pre-acquisition JSON remains the immutable authority.
The purpose is source acquisition, endpoint validation and arithmetic comparison;
**no predictive model fitting is authorized**, even if recovery succeeds.

## Exposure and purpose

The older Zenodo S9 scores have already been fully inspected, modeled and
examined through a separately frozen arithmetic audit. They fail the expected
unweighted ESS identity but satisfy a weighted rank-two constraint. That source
finding motivated this recovery. Old S7 design/reference metadata and the v2
admitted panel are also known. Those results and frozen protocols remain
unchanged.

Before this freeze, only final-attachment metadata, article/preprint methods,
old-source headers and publicly available downstream code were inspected. The
final measurements may be another release of the same experiments. They cannot
be described as an independent cohort, prospective experiment or untouched
predictive holdout. Downloading the complete archive counts as outcome
acquisition, even if only selected members are subsequently interpreted.

The source is the final article, DOI
[10.1038/s41586-024-08430-9](https://doi.org/10.1038/s41586-024-08430-9),
PMCID PMC11903340. Its official MOESM5 caption is “Supplementary Tables 1–13.”
The article declares CC BY 4.0; retain attribution and attachment-specific terms.
The official article XML SHA-256 is
`41a554f31ea67adb1ce983b05ed9194e8fa1596bd4e6d2208df1854de9f33320`.

## Exact acquisition scope

The primary payload is the official
[publisher attachment](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-024-08430-9/MediaObjects/41586_2024_8430_MOESM5_ESM.zip),
`41586_2024_8430_MOESM5_ESM.zip`, exactly **342234089 bytes**, with a 1-GB download
limit. The verified HEAD advertises byte ranges, opaque ETag
`cd38d6f0835534f0227865fee990cb9c` and last modification
`2025-02-14T10:46:05Z`. The ETag is not a promised MD5.

Compute the complete archive SHA-256 after acquisition and before interpretation.
There is no authenticated publisher digest presently available. For range-only
retrieval, report member integrity and never claim the whole archive checksum
was verified. Normal official PMC, Europe PMC supplementary-file and PMC OA
package routes may distribute the same named attachment; preserve package,
redirect, size and checksum receipts. A discrepant copy is a version finding,
not an interchangeable favorable substitute.

Inventory the ZIP central directory before opening values. Validate member
bounds, sizes, CRCs and hashes. Quantitative interpretation is restricted to
the joint-library measurement/specificity fields in these exact members:

```text
2023-03-03628-s5/2023-03-03628C-Table S10 - joint lib design.xlsx
2023-03-03628-s5/2023-03-03628C-Table S11 - jointdata.xlsx
2023-03-03628-s5/2023-03-03628C-Table S12 - folds and performance.xlsx
```

Other table filenames may be inventoried, but their quantitative values are not
part of this recovery. Predictions and performance summaries in any selected
workbook are not empirical activity labels. Save a header/schema-binding artifact
before quantitative checks, naming each worksheet, ID, cell, replicate,
aggregate, ESS, DNA and missing-value role. Missing or ambiguous fields are
nonestimable; do not substitute a surrogate.

## Pre-acquisition code evidence

A downstream MIT repository,
`TaykhoomDalal/LentiMPRA-Enhancer-Processing`, snapshot
`d2751a6ddff67ce4406288e9f6c3d932b0c12cdd`, contains `validate_joint_ess()`.
Its function already existed in the initial public commit dated 2026-08-21.
This predates our recovery but is not producer historical documentation.

The downstream code predicts these bytes before official acquisition:

| Object | Bytes | Predicted SHA-256 |
| --- | ---: | --- |
| Whole archive | 342234089 | `aefb1974e9ebe5c196a89a00f4e868c0abaf69e3d3952e79c74bfb42f217f6a4` |
| S10 design | 8301566 | `f9ec4c6c2d230b0ef507a39e5d5c73b3b36f90ccc5758c33b442c7561e20e7e5` |
| S11 joint data | 53557805 | `9c95ca2b8df0ff3edc7dda4fdcf337f33e652b96a0833a70e786bbcdc2b76b4e` |
| S12 folds | 3564742 | `0f9ff22d4b46fab36d419399e6ab743dc08a0813f812dc7bd96f90eb7e0f2b8a` |

These are **downstream predictions, not publisher-provided checksums**. Their
agreement or disagreement with official bytes is a recorded provenance result.
The predicted S10/S12 identities would make them byte copies of the earlier
Zenodo S7/S9, not independent evidence from new experiments.

The code suggests S11 worksheets `{cell}_replicate_level`,
`{cell}_summary_data`, `all_cell_types_summary` and `all_cell_types_ESS`.
Candidate activity columns are `<cell> [log2(rna/dna)]`; candidate ESS columns
are `ElementSpecificityScore_<cell>`. Verify these against the actual official
headers rather than silently assuming their meaning. The downstream project
explicitly calls the upstream normalization of S11 ESS unpublished.

The separately identified GAME evaluator repository remains **metadata-only**
under this protocol: it has no declared repository/data license and a missing
intermediate file. Its workbook payloads are not authorized sources here.

## Endpoint and row accounting

Keep all original empirical values, source numeric representations and missing
tokens unchanged. Bind observed joint-library log2 RNA/DNA activities and their
documented three infection replicates separately from published aggregates,
ESS, counts and predictions. Replicates belong to one source study/library;
barcodes, molecules and repeated file copies are not independent elements.

The published methods describe per-replicate median normalization followed by
averaging three replicates. They describe per-cell activity z-scoring followed
by per-element across-cell centering for ESS. These are semantic claims to
verify. Do not infer unavailable raw activities, pseudocounts, replicate values
or normalization constants from ESS.

Use exact case-sensitive stable IDs for primary old/new comparisons. A global
documented naming alias may be bound only through an ID-only audit before
quantitative comparison; retain exact-join results separately. No fuzzy matching
or choice of mapping by score agreement is allowed. Preserve source row IDs,
duplicates, blank IDs and missingness. Conflicting duplicate observations block
the affected comparison rather than being averaged away.

Report all source IDs and each maximal complete three-context endpoint cohort,
then their exact old-S9 intersections and intersections with the fixed 55276
previously admitted IDs. No activity, ESS, p-value, category, count, rank or
agreement filter is added. Missing required replicate values make that check
nonestimable; no `nanmean` or imputation is allowed. Normalization populations
are fixed before these comparison intersections, never optimized for agreement.
Compare physical DNA/coordinate identity when final fields are available.

## Finite arithmetic checks

All derived arrays are diagnostics, never replacement targets. If their fields
exist, compute these three activity bases:

1. The unchanged documented published aggregate.
2. The arithmetic mean of three finite published replicate activity scores.
3. The arithmetic mean after subtracting each replicate's median over all its
   finite source values. Record every median and population count.

For each base, report the complete fixed grid:

- Pre-z-score population: all complete three-context source rows, or each
  context's own finite source population followed by complete-row evaluation.
- Pre-z-score SD convention: population `ddof=0` and sample `ddof=1`.
- Subtract each element's three-context mean.
- Post-step: none, column z-score with `ddof=0`, or column z-score with `ddof=1`,
  using the full complete source population before old/admitted intersections.

This gives **36** potential variants. Report every one, including equivalences
and nonestimable cases. The extra column-standardization variants are explicitly
forensic compatibility checks motivated by the already observed S9 structure;
they are not established historical steps. No continuous transform search,
fitted affine calibration, optimized weights or favorable subset is allowed.

There is one additional, separately specified pre-value code hypothesis:

1. S11 published `all_cell_types_ESS` obeys the unweighted zero-sum identity.
2. For every exact S12 score ID, S12 observed ESS equals the corresponding S11
   ESS minus its mean on the **entire exact S12 ID subset**, divided by sample
   SD on that same subset (`ddof=1`).

Require every S12 ID to have a unique finite three-context S11 match. Record
the ID-set hash; missing matches make the test nonestimable rather than defining
a smaller favorable subset. The maximum absolute tolerance is **`1e-12`**.
Report the full S12 cohort and the fixed old-S9/admitted intersections. This
checks a release relationship, not producer history.

Compare each estimable replicate aggregate to the published aggregate, final
published ESS directly to unchanged old S9, and every fixed pipeline to final
ESS and old S9 where present. Keep exact Decimal equality counts separate from
floating agreement. The general floating rule is

`abs(candidate-source) <= 1e-8 + 1e-8*abs(source)`.

Report matched/missing counts, signed mean error, MAE, RMSE, maximum absolute
error, Type-7 absolute-error quantiles, all context means/SDs and ESS row sums.
Do not widen tolerances or omit unfavorable variants.

## Interpretation and completion

Maintain separate flags for source authenticity, empirical endpoint binding,
replicate arithmetic, ESS formula, old/new release relationship and historical
pipeline identification. A numerical match establishes compatibility, not
necessarily the unique historical process. Positive weighted geometry alone
cannot identify normalization; per-cell offsets/scales can disappear under
z-scoring, and some SD conventions become indistinguishable after a subsequent
column standardization. Preserve such nonidentifiability explicitly.

Blocked, missing, divergent and unresolved results are valid recovery outcomes.
No old target, protocol, fit or reported loss is changed. Any future predictive
analysis needs a verified endpoint contract and a separate frozen analysis
protocol; this recovery does not authorize it.
