# Testing whether activity-model slopes transfer to variant effects

## Result: the frozen v1 source gate failed

Both source readers stopped at the same first failing selected effect cell:
**`AO42972`**, K562, variant **`chr12:72652619:A:G`**. Its literal value is the
shared-string text **`NA`**. The committed v1 binding accepts no missing-value
token, so the source observation gate failed. This closes v1 with a
**source-validation failure and no model-accuracy result**.

The independent gate comparison matches all **241 selected-cell ledger rows**,
including exact binary64 bits for **240 ordinary numeric cells** representing
**120 complete paired variants** before the stop. That agreement reproduces a
failed gate; it cannot turn it into a pass. No predictive-performance metric,
bootstrap interval or favorable partial-pair analysis was computed.

The pre-outcome metadata audit had fixed **2,172 human SNVs in 172 components**
for K562 and HepG2. Its independent grouping verification remains valid:
**4,689 comparison categories** cover all 1,171 upstream components and every
one of the 332,719 reference-validated variant assignments, prior-exposure flags
and exclusion reason classes. The metadata cohort never became a verified
panel of 4,344 finite reporter observations.

The primary attempt stopped before making predictions. The independent attempt
had already exported fixed DNA/model predictions before it reached the same
source failure; those predictions are retained, without a measured accuracy
claim. The [failed-v1 evidence snapshot](results/variant-transfer-v1/README.md)
preserves both attempts, the selected-cell ledgers and their comparison.

This distinction matters because the [earlier activity result](context-activity-study.md)
concerns held-out **whole-fragment activity** within the Agarwal experiment.
Its 8.7996% MSE improvement does not establish transfer to single-base effects,
a separate assay or an independent biological mechanism.

## Frozen question and exposure history

The authority is [variant-transfer-v1](../experiments/variant-transfer-v1.json),
SHA-256 `811930c2d3a5ee763014ea611acda4d2c84404b7d5ae399c7843c62d0a828149`;
[the readable protocol](../experiments/variant-transfer-v1.md) accompanies it.
The semantic contract was frozen in commit
`ec454e6d559dc96337359e32f8a05242581798c2` **before the first S4 acquisition
request**. The [official acquisition receipt](variant-transfer/source/siraj2026-S4-acquisition.json)
records that request and its protocol binding. This is a repository freeze,
not external preregistration.

Earlier HBB/Kircher outcomes, Agarwal specificity/activity values and fitted
models, source-selection discussions, and the new article's methods and
aggregate metadata were already known. DNA-only S2 had been acquired and its
schema inspected before the semantic freeze. A later DNA check corrected the
nominal 336,475 primary-chromosome SNV-pair count to **336,471 literal ACGT
pairs**: four nominal pairs contained ambiguous bases. No reporter outcome
selected that correction.

The design separates two freezes. The first fixes source, endpoint semantics,
contexts, population, predictors, losses and uncertainty before downloading S4.
The second binds exact textual headers, source preprocessing, missing policy,
DNA identities, all prior-exposure/control exclusions, original global component
IDs and fixed model hashes before parsing numeric effects. That second commit is
`a21fd6f44f493b20cddc51e6cde7c7ae14e490e4`; it binds
[the final metadata receipt](variant-transfer/binding-v1-verified/binding.json),
SHA-256 `9c26b3fc628f16bd6bb677dee8eed01e509b6102d44f95c0037ae2f55676ba79`.
Its `S4_numeric_effect_values_read: false` describes the exposure state at the
freeze, before the subsequent failed observation attempts.

Metadata-only inspection did not authorize choosing a different target or
population. The numeric failure did not alter either frozen rule set. A semantic
change requires a new exposed protocol version.

## Source, endpoint and assay limits

The external experiment is Siraj et al., *Functional dissection of complex
trait variants at single-nucleotide resolution*, Nature (2026), DOI
[10.1038/s41586-026-10121-6](https://doi.org/10.1038/s41586-026-10121-6).
The [source packet](variant-transfer/source/manifest.json) retains official
receipts, methods, design/header checks, producer-code metadata, assembly
mapping and the source-semantics review.

| Asset | Role | Bytes | SHA-256 |
| --- | --- | ---: | --- |
| Official S2, `41586_2026_10121_MOESM4_ESM.xlsx` | Exact inserted oligo DNA | 40,795,798 | `4b6e20d1789021ee8f11f57af366a65680d25496de9351a9e781e9424716d824` |
| Official S4, `41586_2026_10121_MOESM6_ESM.xlsx` | Published main-panel results | 224,405,074 | `723c4a0e1cfdb2c9c6ae55b3fa9206855dd95e0241319b9fc8d0831d3b42f9c2` |

S4's `S4.MPRA results` sheet has header row 2, ID column A `Variant`,
AN `log2Skew in HEPG2`, and AO `log2Skew in K562`. The
[header-before-values receipt](variant-transfer/source/S4-header-inspection-v1/headers-before-values.json)
and [source-semantics review](variant-transfer/source/source-semantics-review.json)
bind the endpoint as the main-panel **DESeq2 Alt-versus-Ref by RNA-versus-DNA
interaction**, in published summit-shift-normalized delta log2 fold-change
units. Positive means greater RNA-relative-to-DNA activity for the alternative
construct. Sign is established from methods and allele identity, not selected
from observed performance.

These are source-estimated allelic effects in a plasmid reporter with a minimal
TATA/GFP design. They are not absolute transcription, endogenous target-gene
expression, a raw-count ratio or a literal subtraction of the earlier Agarwal
aggregate activity estimates. The source saturation-mutagenesis posterior
endpoint has different processing/sign concerns and is excluded; S7 multivariate
shrinkage, source prediction columns and across-cell meta-analysis are also
excluded. The two named cell lines do not make assay delivery, culture or
normalization identical between studies.

The source methods exclude variants with fewer than 30 DNA counts or exactly
zero RNA counts for either allele. For repeated libraries, they choose by
plasmid coverage; S4 reports the source-selected measurement per variant/cell.
The test uses that authoritative published preprocessing unchanged. It does not
reconstruct discarded libraries, apply guessed thresholds to rounded aggregate
counts, or filter by activity, emVar, FDR, PIP, CRE annotation or effect size.
ID uniqueness and methods provenance do not independently audit the raw counts.

Publisher article licensing is recorded as CC BY 4.0, with source-workbook
terms reviewed separately. The linked Zenodo deposit's CC BY-NC 4.0 terms are
not silently substituted for the official publisher distribution. Large S4
bytes remain in a checksum-bound cache outside Git; receipts and exact selected
bindings are retained in the repository.

## Outcome-free population, assembly and exclusions

Both source alleles must contain an exact **200 nt uppercase ACGT genomic
core**, differ at exactly one position and agree with the source ID's Ref/Alt
bases at zero-based index 99. No indel, inferred/padded core, ambiguous base,
synthetic sequence or outcome-selected active subset is admitted. Documented
adaptors are not treated as genomic DNA.

The source variant positions are in hg19. Each full insert is mapped through a
checksum-bound hg19-to-hg38 chain within one aligned block and then compared
exactly with the locked hg38 reference, including strand. Merely lifting a
single variant coordinate is insufficient. Missing/ambiguous mappings, primary
chromosome failures and full-reference mismatches are retained in the design
ledger. No hg19 coordinate is relabeled as hg38.

| Metadata stage | Variants | Components |
| --- | ---: | ---: |
| All unique S2 design IDs | 365,762 | Not yet formed |
| Literal primary-chromosome 200 nt SNV pairs | 336,471 | Not yet formed |
| Full-insert chain mapping and exact primary hg38 reference match | 332,719 | Full population used for grouping |
| Complete global sequence/locus partition | 332,719 | 1,171 |
| After whole-component prior-exposure exclusion | 2,430 | 174 |
| Final source-ID/control metadata cohort | **2,172** | **172** |

The [design mapping summary](variant-transfer/source/design-mapping-v2/summary.json)
retains 1,040 pairs without a whole-insert single-block mapping, 2,700 full
reference mismatches and 12 mappings to nonprimary hg38 contigs. Its earlier
rowwise `eligible-pairs` file is a diagnostic, not the final grouping population.
All **332,719** reference-validated pairs enter the global graph before any
old-source exclusion.

That graph joins source parents, exact/RC allele identities, positive genomic
interval overlaps, **floor(start/1,000,000)** bins keyed by chromosome, and full
200 nt ungapped at-least-90% identity between either allele on either strand.
Links are transitive. Crossing a bin boundary does not itself join every touched
bin; actual interval overlap is a separate relation. All contexts of a variant
retain the same original component. No component is reconstructed from a
prefiltered set or regrouped after exclusions.

The [prior-exposure inventory](variant-transfer/prior-exposure-inventory/README.md)
retains all 55,276 Agarwal training cores, all fifteen explored Kircher assays
at ten promoter loci, their source-resolved substitutions/deletions, and HBB
reference/ClinVar intervals. A whole new component is excluded if any member
lies within the inclusive 10,000 bp old-locus guard or has qualifying old-DNA
homology. Old nodes never receive bins that could merge otherwise separate new
components. All old alternate-contig cores remain in DNA screening, while new
Siraj records remain restricted to primary chromosomes.

Variable-length old promoter targets are retained in full. Their derived
200 nt reference windows, factorized observed SNV alternatives and explicit
deletion-created windows support exhaustive old-DNA screening. The isolated
edits are benchmark molecular identities, not proof of original single-mutant
barcode constructs. HBB's shorter promoter remains covered by its genomic
locus guard. No equality-form HGVS record is called neutral: original author
alleles resolve the previously explored deletion identities.

The official S30 control design supplies **41 exact/RC control-matched variant
IDs** in **three global components**. Each affected component is excluded.
Those components were already prior-exposed, so this does not remove an
additional member of the 2,430-member prior-cleared set. Source ID-pattern
concordance is a diagnostic; sequence identity supplies the control relation.
Genuine newly measured location/annotation/null controls remain eligible.

S4's ID-only inventory has 304,278 unique IDs, all present in S2. Of the
332,719 reference-validated pairs, 54,468 are absent from S4. Intersecting the
prior-cleared population with authoritative S4 membership removes **258**
variants and leaves two components without an admitted member, producing the
2,172/172 metadata cohort. Its original component IDs are preserved. This
membership check reads IDs, not activity, allelic effects or significance.
Exclusion reason counts overlap and must not be summed as disjoint categories.

## Independent cohort reconstruction

The primary Python grouper and an independently written C++/Python
implementation start from the full mapped DNA and old inventory. They use
different packed-DNA representations, index data structures and genomic-guard
algorithms. Twenty-one disjoint exact pieces guarantee candidates within twenty
mismatches; the factorized old SNV branch uses twenty-two pieces to retain the
twenty-one-mismatch cases that a documented old substitution can repair.
Every candidate is checked exactly. Pruning removes only redundant unions or
further old-hit searches for an already excluded whole component.

The [independent comparison](variant-transfer/independent-grouping-v1/comparison-v1/comparison.json)
passes **4,689 categories with zero failures**, including all component members,
sizes, every variant's retained/excluded status and reason classes, and both
input-source bindings. All 1,656 homology union-witness triples also match.
Both complete output manifests are authenticated before comparison. First
old-hit traversal choices and coordinate-witness representation are retained
in [witness history](variant-transfer/independent-grouping-v1/witness-history.json),
without requiring an arbitrary identical first witness.

Fourteen executable independent controls cover bin/overlap boundaries,
chromosome namespaces, source parents, Hamming thresholds, reverse strands,
SNV repair, inclusive guards, whole-component propagation and old nodes
remaining separate. This verifies software reconstruction of the metadata
cohort; it does not supply an independent biological observation or validate
numeric target coverage or measure predictive accuracy.

## Fixed predictor, loss and uncertainty

The candidate is the **matched-cell I model**, averaged across the five fixed
Agarwal activity fits. The clean activity core's artifact-manifest SHA-256 is
`2ea476bc5413dc13486dec9f3d7e9162a50a0e47f1d83ce29efabc5f5cc5af9d`,
from code commit `035ca5731fe594622cd2c3c27cd3bb51ec65c924`. No new data fit,
full-training refit, response calibration, sign selection, coefficient editing
or penalty tuning is authorized.

Features remain the existing 32 RC-canonical overlapping 3mer frequencies plus
GC fraction, with denominators 198 and 200. For each cell:

```text
predicted_delta = mean_five_fits sum_features (
    standardized_coefficient / old_training_feature_scale
    * (feature_alt - feature_ref)
)
```

Intercepts and old training-feature means cancel. The baseline predicts **zero
effect**. Shared A/S slopes provide a fixed secondary predictor; their pair
differences coincide. There is no newly trained GC-only baseline. GC-preserving
versus GC-changing substitutions are metadata-defined diagnostics.

For a component `g` with admitted variants `i`, define:

```text
L_g(model) = mean_variants_in_g mean_{K562,HepG2}
             (observed_delta - predicted_delta_model)^2
primary_gain = mean_components [L_g(zero) - L_g(I)]
```

Each of the two cells receives equal weight within variant; each admitted
variant receives equal weight within component; each component receives equal
primary weight. The number of admitted members, rather than the original
pre-exclusion component size, determines its within-component mean. Positive
macro-component MSE improvement favors the fixed model. Evidence of useful
transfer requires positive observed improvement and a positive lower endpoint
of the declared conditional interval. Secondary per-cell losses, MAE, equal-
variant sensitivity, signed correlation or sign agreement cannot replace a
failed primary comparison.

The unexecuted uncertainty plan specifies **2,000 paired whole-component draws**,
sampling with
replacement from lexicographically ordered component IDs using verified
SplitMix64 seed `0x2026102801`. Each draw carries all admitted variants and both
contexts, with fixed predictions. It calls for Type7 2.5%/97.5% quantiles and all
draw indices and losses. Those would be conditional selected-panel ranges, omitting
old-model refitting, source normalization uncertainty, within-assay estimation
covariance and study-to-study variation. They are not calibrated population
confidence guarantees or evidence of a causal regulatory mechanism.

## Independently reproduced observation-gate failure

The [verified metadata binding](variant-transfer/binding-v1-verified/binding.json)
fixes `missing_policy.allowed_tokens = []`. The
[missing-token review](variant-transfer/source/missing-token-evidence.json)
found no authoritative literal missing-value convention for the final S4
workbook. Generic R-writer behavior was not substituted for missing export
provenance. Required blanks, absent cells, formulas, errors, booleans, numeric
strings, unrecognized text and nonfinite values therefore block this v1.

After the second commit, the primary strict-XML reader and the independent
openpyxl reader both reached `AO42972` for `chr12:72652619:A:G`, K562. The XML
cell type is `s`, and the decoded literal text is `NA`. Neither reader silently
converted it into a documented missing observation, dropped the variant,
imputed a value or continued to a smaller favorable cohort. The observation
gate is **false** under the unchanged v1 policy.

The comparison receipt originally written at
`runs/variant-transfer-v1-gate-comparison/comparison.json`, retained in the
[failed-v1 snapshot](results/variant-transfer-v1/README.md), has SHA-256
`a13b8ea6aaf77a473ee50898e0d35031324af5bba6932a9a554d5ab013fc38fc`.
It authenticates both attempt manifests and reproduces every one of the 241
selected-cell ledger rows through the failure, with **zero comparison failures**.
All 240 ordinary numeric source float encodings agree exactly. These form
120 complete pairs before the unsupported selected cell; they were not used as
a substitute accuracy-evaluation population. This is a partial source-gate
comparison, not a full successful protocol replay.

The primary attempt produced no predictions. The independent reader exported
`fixed-predictions.json.gz` (**353,893 bytes**, SHA-256
`b6552fd7e5fa633c65188f8faeefe31c4105e211dec2113c184fd8a4ae27ceeb`)
from the frozen models and DNA before source parsing failed. That artifact is
preserved as part of the attempt history. Its existence is compatible with the
absence of any measured accuracy result: **no MSE, MAE, correlation, sign
accuracy, bootstrap or transfer-performance claim was computed**.

The readers may differ in internal parsing: openpyxl can instantiate
nonselected worksheet values while iterating. The ledger comparison audits
selected-cell use through the first failure, not a guarantee that no other
cell value was ever represented internally. Nothing in this software agreement
establishes the biological meaning of the unsupported string or verifies full
numeric coverage.

At least 100 complete variants and ten components were required for an
inferential claim. The metadata counts exceeded those minima, but the failed
source contract prevents invoking that support threshold or evaluating the
remaining performance rules. V1 remains a source-validation failure rather
than evidence for or against the model's accuracy.

Two grouping implementation attempts remain visible:

1. The first preparation mistakenly encoded all touched 1 Mb bins. The frozen
   protocol instead specifies start bins plus positive overlaps. Its
   [error receipt](variant-transfer/source/grouping-preparation-error-v1.json)
   preserves the preparation and records that no full run or numeric outcome
   parse occurred. The implementation was corrected to the existing contract;
   the protocol was not amended to match the mistake.
2. The first corrected full-run attempt applied a primary-chromosome validator
   to old prior intervals and stopped before grouping. Four valid old training
   cores lie on `chr22_KI270879v1_alt` or `chr7_KI270803v1_alt`. The
   [failed run](variant-transfer/source/grouping-v1/failure.json) is retained.
   The adapter was corrected to retain every old core and its DNA while keeping
   new variants primary-only. The successful
   [primary grouping](variant-transfer/source/grouping-v1-prior-contigs/summary.json)
   and independent reconstruction use the unchanged full old inventory.

Neither attempt inspected reporter effects or selected a source transform,
cohort rule or tolerance to improve predictive performance. Earlier metadata
binding preparations also remain alongside the final verified binding; they
are not separate fitted experiments.

## Separate exploratory continuation and scope

An operational numeric-complete-case continuation belongs to the separately
versioned [exploratory v2 protocol](../experiments/variant-transfer-v2.json),
with its own freeze governing further parsing and analysis. It follows exposure
to the v1 failure, 120 complete observed pairs and precomputed fixed predictions.
It cannot be described as an outcome-naive continuation or a retroactive v1
success. No v2 performance result is reported here.

V1 supplies an inspectable cohort/source audit and a reproducible stopping
result. It does not establish transfer accuracy for this small composition
model, a globally new modeling method, whole-genome understanding, endogenous
target-gene causality, clinical consequence or an independently established
biological mechanism. Every earlier v1 rule, source failure and implementation
attempt remains retained.
