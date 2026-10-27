# Original published joint-library source recovery

This directory records acquisition, integrity, exact member selection, a schema
gate written before values were decoded, and an unchanged source-value export.
The governing protocol is `experiments/context-endpoint-recovery-v1.json`, SHA256
`00175e41cf9fa0cbf177ce8fc7ef143eb508aec8b11d058f7acb4458284dc868`. Parent committed
that protocol before the new attachment GET. No predictive fitting or activity
arithmetic is performed by the acquisition/export tools here.

The official final Nature supplementary ZIP was obtained from the exact
publisher link recorded by the metadata discovery round. It contains 13 members,
has 342,234,089 bytes and SHA256
`aefb1974e9ebe5c196a89a00f4e868c0abaf69e3d3952e79c74bfb42f217f6a4`.
The complete download receipt is `official-archive-acquisition.json`; parent
validated the CRC of every member in `archive-validation/verification.json`.
There is no publisher-provided digest: the observed SHA256 matches a downstream
code prediction recorded before download, which is a separate fact from official
publisher provenance. The HTTP ETag is not asserted to be MD5.

`archive-inventory.json` lists every literal archive member, compressed and
expanded sizes, local offsets and CRCs. Only the three exact joint-library
members permitted by the protocol were extracted. `extract_members.py` checks
protocol/archive hashes, selects unique literal names, bounds expanded bytes,
rejects traversal, symlinks and encryption, checks CRC at EOF and publishes
without overwriting an existing file. Its synthetic tests include a malicious
pre-existing output ancestor symlink and prove that no directory is created
outside the chosen output tree before rejection.

`member-extraction.json` and `member-source-identities.json` establish:

| Final article table | Source workbook | Retention and identity |
| --- | --- | --- |
| 10 | joint library design | Byte for byte identical to the older retained Zenodo S7 workbook |
| 11 | jointdata | Retained unchanged at `datasets/context-endpoint-recovery/raw/jointdata.xlsx`; 53,557,805 bytes; SHA256 `9c95ca2b8df0ff3edc7dda4fdcf337f33e652b96a0833a70e786bbcdc2b76b4e` |
| 12 | folds and performance | Byte for byte identical to the older retained Zenodo S9 workbook; its observed-score sheet is reused, and benchmarking outcomes are not part of this recovery |

The large outer ZIP stays in the runtime cache. The existing retained S7/S9 files
are reused rather than duplicated. Exact byte identity includes all original
design information and observed-score representations; the copies are not
independent measurements. Different preprint/final table numbers do not imply
different experiments.

`headers-before-values.json` inspected first-row **text headers only**. Numeric
first-row cells were redacted. `replicate-label-metadata.json` inspected only the
replicate-label column and established literal labels `1`, `2` and `3` before
binding. `schema-binding.json` binds all eight worksheets and writes explicit
machine-readable aggregate, long-form replicate, published ESS and separate
summary-replicate roles. Its SHA256 is
`a59592816f57f7254b780c83ba3e929e3b38f4fd06b2b196a0cddda1d9bfb98e`.
The HepG2 replicate sheet has additional numeric cells in first-row columns
H–K; these remain unbound and uninterpreted. No activity or ESS values were
decoded before the binding was saved.

The primary published activity aggregate is the explicitly labeled
`all_cell_types_summary` sheet. The primary replicate activity uses the
`log2(rna/dna)` field of each cell's `replicate_level` sheet and its explicit
replicate label. The per-cell `summary_data` replicate fields and `mean` are
retained separately. Selection was based on these source roles before numerical
checks; no rounded/unrounded choice was made by agreement. Headers do not prove
which median normalization has already been applied to each field. That question
belongs to the frozen arithmetic audit, as does the relationship between the
separately labeled `all_cell_types_ESS` sheet and the previously retained S9.

`export_bound_values.py` verifies all source/binding/protocol hashes and writes
`jointdata-source-values.json.gz`. Tables use exact source sheet names and
headers, with each row retaining `id`, the one-based XML `source_row` and
`values` containing unchanged numeric lexemes or original missing tokens. The
export preserves lexical precision; it does not calculate, round or replace
observations. `source-type-missing-formula-ledger.jsonl.gz` preserves non-default
cell types, absent/missing/nonfinite/text cells, cached formulas and their
attributes, plus duplicate/missing ID and empirical-key ledgers. Default finite
numeric cell types are counted per field; their complete original lexemes remain
in the compact export. `source-export-receipt.json` binds output hashes and
counts. Repeated element IDs in long-form replicate tables are distinguished
from repeated `(element ID, replicate)` empirical keys.

The shared XML row decoder is reused from `scripts/fetch_context_panel.py`.
`openxml_source.py` adds explicit ZIP/shared-string bounds and duplicate-member
checks while streaming worksheet rows. Both implementation hashes are recorded.
The actual workbook expands to 265,800,125 bytes. The raw XLSX remains the
authoritative source for omitted count, ratio, barcode, category and formatting
fields. Omitting those fields from the compact activity export does not create
a new source-QC filter.

Source: Agarwal et al., *Massively parallel characterization of transcriptional
regulatory elements*, Nature (2025), DOI
[10.1038/s41586-024-08430-9](https://doi.org/10.1038/s41586-024-08430-9).
The final article grants CC-BY-4.0, including supplementary material except where
a specific third-party credit line says otherwise. Its exact permissions and
official links are retained in `research/context-provenance-round3/`. Preserve
attribution, license and transformations. No conflicting terms were identified
in the inspected source metadata/headers; arbitrary data-cell annotations were
not searched for licensing assertions. Receipts retain safe response headers
only, without cookies, nonces or opaque reporting fields.

Earlier v1's missing raw-activity semantic gate and all earlier frozen results
remain unchanged. Acquiring this final source does not repair an earlier
protocol retrospectively. Three infections belong to one source study/library;
these files do not create an independent validation cohort.
