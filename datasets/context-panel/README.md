# Reference-validated three-context specificity panel

This dataset contains **55,276 deposited DNA elements × three published Element
Specificity Scores (ESS)** from the preprint-associated Agarwal joint-library
release at [Zenodo10558183](https://zenodo.org/records/10558183), related to
[Agarwal et al., Nature2025](https://doi.org/10.1038/s41586-024-08430-9).
It contains a selected source panel, not a random sample of human regulatory DNA.

**Source-formula validation remains unresolved.** The frozen workbook values
do not satisfy the unweighted across-cell centering formula described by the
papers. Exact import/reference validation and numerical model reproduction do
not resolve that discrepancy. See the
[semantic review](../../research/context-specificity-semantic-review.md).

## Two distinct protocols and an explicit failed acquisition

`context-prediction-v1` required measured mean log2 RNA/DNA activity. The retrieved
S7 workbook contained design metadata only; S9 explicitly contained ESS, which
v1 forbids as a raw-activity label. **V1 stopped without fitting.** The failure,
source schema and design-only attempt remain in `v1-*.json` and the losslessly
repacked `design-catalog-v1.json.gz`. Its original and compressed byte commitments
are in `v1-design-artifact-repacking.json`.

The separate frozen `context-specificity-v2` protocol, SHA-256
`8f0576845fecb5cfcf9bfa354fab446d849244448a3489905a6ce33bdeb21cec`,
authorizes evaluating the published ESS instead. The workbooks were already
retrieved, but only header/ID/design metadata had been inspected before v2 was
frozen. The importer records this acquisition history and never relabels the v1
failure as a successful activity analysis. Source CNN prediction files were not
used as measurements.

## Endpoint and interpretation

The papers describe z-scoring log2 RNA/DNA activity within each cell, then
subtracting each element's mean across HepG2, K562 and WTC11. That formula
requires the three resulting scores to sum to zero. The frozen S9 columns
violate it: the maximum absolute sum is about1.09873, far above rounding error.
They instead satisfy a precise weighted rank-two relation with unequal cell
weights. Extra column scaling is a compatible explanation, not established
producer history. These dimensionless published score columns are retained
unchanged; they are not raw-expression or log2-activity units.

The documented target processing uses full-panel statistics, including future
held-out elements; the workbook's historical transform is unresolved. No
prospective or training-only target preprocessing is claimed. Under the stated
formula, zero specificity does not mean inactivity and common activity is
removed. Those interpretation rules cannot silently establish S9's historical
normalization. Values are not recentered, repaired or forced to sum to zero.

The original `source-schema-binding-v2.json` acceptance flag records the exact
worksheet/column-role gate before this arithmetic discrepancy was discovered.
It is not full formula validation. Its bytes remain preserved; the later
semantic review and published run explicitly retain semantic validation=false.

The protocol controls downstream features, grouping, fits and losses. This
importer does not fit models or select records by score magnitude, predictability,
source category, rank or agreement with a model.

## Source and representation

The untouched source workbooks in `raw/` retain their original published byte
sizes and verified Zenodo MD5 checksums, plus observed SHA-256 hashes:

| Workbook | Role | Size bytes |
| --- | --- | ---: |
| `Table-S7-joint-design.xlsx` | 60,000 designs, genomic coordinates and230nt constructs | 8,301,566 |
| `Table-S9-folds-performance.xlsx` | Published per-cell ESS with source IDs | 3,564,742 |

The source title says S7/S9. They are not silently renamed as final-paper
Supplementary Table10 or asserted to be byte-identical to final publication files.
Zenodo declares CC-BY4.0 for its record and listed files; the exact metadata grant
is retained as `raw/zenodo-record-metadata.json`. Preserve the source attribution.
A code repository's software license is not substituted for this data grant.

S7 explicitly labels the230nt construct as15nt adapter +200nt genomic element
+15nt adapter. The core is exactly `sequence[15:215]`; no motif or reference
search chooses another subsequence. Full constructs, adapters, all component
hashes, strand and source-row provenance are retained in
`source-design-metadata-v2.json.gz`. The normalized panel uses the200nt core.

S9 has55,338 nonempty unique source IDs, all exact unique S7 matches. No `_F`
suffix, prefix, fuzzy match or outcome-dependent ID transformation is used.
Its worksheet dimension ends at row59108, but only55340 physical row elements
exist. A maximum row number is not a biological record count. The blank-row/gap
and ID-only audit is retained rather than manufacturing missing element IDs.

## Independent reference and complete inclusion accounting

Every admitted core exactly matches the independently downloaded UCSC hg38
reference at its **declared** zero-based half-open interval and source strand:
54,195 forward matches and 1,081 reverse-complement matches. The reference SHA-256
is `1f67aaa17a77b327738fe750ab37430a85dddf73c7d9189385ad1839259acec0`.
Declared verified alternative contigs are eligible; there is no unstated
primary-chromosome-only restriction, silent coordinate search, liftover rescue,
strand correction or sequence replacement.

The importer accounts for all source designs and IDs:

- 55,276 complete, finite, exact-reference elements are admitted.
- 52 exact-joined S9 IDs have reference-core mismatches;10 have non200nt intervals.
- One explicit S9 row has no stable ID.
- 4,662 S7 designs have no exact S9 ID, including source controls or other designs
  absent from this published score table. They are not silently imputed.
- The complete60,000-design audit separately records670 missing-coordinate
  designs,11 non200nt intervals and 55 reference mismatches; these larger counts
  include designs not present in S9.

`inclusion-report.json` retains all source exclusions, missingness and per-context
numeric status accounting. Raw coefficients/ESS are never dropped merely because
their values are inconvenient. Distinct source IDs are preserved even for
identical or reverse-complement sequences. No undocumented biological parent
relationship is invented; ROOT-owned grouping subsequently combines exact DNA,
parent/overlap/bin and declared homology relationships under the frozen protocol.

## Artifacts and verification

| Artifact | Purpose |
| --- | --- |
| `panel.json` | Exact admitted core DNA, coordinates, source IDs and three ESS observations |
| `manifest.json` | Protocol/endpoint/source/normalization commitments and every consumed artifact hash |
| `source-schema-binding-v2.json` | Exact worksheets, headers, semantic role and version limitations |
| `inclusion-report.json` | Complete source/missingness/exclusion accounting |
| `reference-validation.json` | Admitted and full-design reference checks |
| `source-design-metadata-v2.json.gz` | Full deposited constructs, adapters, hashes and source-cell provenance |
| `raw/retrieval.json` | Original workbook retrieval/size/MD5/SHA receipts under v1 |
| `v2-ess-inspection-gate.json` | Distinct frozen v2 authorization before numeric ESS parsing |

Integrity, exact source-cell values and design provenance can be replayed offline:

```bash
python3 scripts/fetch_context_panel.py verify datasets/context-panel
```

To replay every admitted interval/strand against the independent cached reference:

```bash
python3 scripts/fetch_context_panel.py verify datasets/context-panel \
  --reference /workspace/.tools/reference/hg38.2bit
```

The output reports whether the full reference was actually replayed. Missing
reference data are not described as verified by an unexecuted check. The script's
`self-test` exercises OpenXML strings, numeric representation, source errors and
header-only column selection without requiring a model or network request.
