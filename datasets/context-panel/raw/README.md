# Frozen post-protocol acquisition receipts

The two untouched workbooks were downloaded only after protocol
`context-prediction-v1` was frozen at SHA-256
`a089a6668ad73e943797df405b3ba510b1377190b416da1bd8df12bc368b77c2`.
Published Zenodo byte sizes and full-file MD5 checksums were verified before
schema inspection. `retrieval.json` retains URL, retrieval time, byte count,
MD5 and SHA-256.

Source: Agarwal et al., *Massively parallel characterization of transcriptional
regulatory elements*, DOI10.1038/s41586-024-08430-9, Zenodo10558183.
The Zenodo data record declares CC BY4.0. Source names and attribution remain
unchanged in the receipts; local filenames are filesystem conveniences.

The empirical-target gate **failed**:

- S7 provides60,000 design records with exact230nt assay sequences, a declared
  200nt core plus15nt adaptors each side, hg19/hg38 coordinates and strand.
- S9 provides **55,338 actual nonempty, unique sequence IDs**, all exact matches
  to unique S7 design names. Its three label columns explicitly contain
  **Element Specificity Score**. Its second sheet contains model benchmarks.
- Neither acquired source contains the protocol's required observed mean-of3
  per-replicate-median-normalized log2(RNA/DNA) label.

ESS, model predictions and per-cell zscores were expressly prohibited as
substitute targets. No fit, performance calculation or endpoint switch occurred.
`schema-inspection.json` retains header inspection with numeric cells redacted.
Its worksheet dimension is **A1:E59108**, not a record count. The sheet stores
**55,340 XML row elements**: one header, 55,338 nonempty unique ID rows and one
explicit blank-ID row. There are 55,339 stored rows after the header.

The original semantic receipt and earlier summary incorrectly called
dimension-minus-header (59,107) the number of IDs. That claim is superseded by
`source-accounting-correction.json`, based on the existing ID-only audit.
`semantic-gate-corrected-metadata.json` records the corrected count while keeping
v1's endpoint failure. Original `semantic-gate.json` and `schema-inspection.json`
remain byte-exact historical evidence; the original README claim is retained
in `README-v1-count-claim-superseded.md`.

Workbook checksums, original acquisition times, source license and v1 protocol
binding are unchanged. This correction read metadata and ID strings only;
no ESS numeric values were parsed or inspected.

Do not treat the presence of these files as a validated model-ready panel.
