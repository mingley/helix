# Activity binding adversarial review

The fixed source/cohort checks pass. One material gate gap was reported to parent
before building: the current builder accepts a recovery summary without checking
that its recorded arithmetic was run on the same export, schema, raw source and
ledger used for activity binding. A stale or incomplete summary from the same
recovery protocol can satisfy the current flag-only test. Parent owns any fix;
this review did not change the builder or tests.

Reviewed code: `scripts/build_context_activity.py`. Its exact review-time hash,
source/cohort checks and license metadata are retained in
`activity-binding-source-checks.json`. The eight existing synthetic controls in
`research/independent-context-activity-replication/test_activity_binding.py`
passed unchanged.

## Requested gate correction

The builder currently checks recovery protocol SHA, `variant_count == 36`, two
status flags and a nonempty aggregate-comparison object. It does not verify a
hash relationship from that summary to the current source export, binding,
raw source, ledger or export receipt. It also does not verify that the claimed
36 variant records are present and unique. This matters because recording the
complete frozen recovery audit is a prerequisite to activity binding and later
fitting; a protocol name alone does not establish which inputs were audited.

The recovery producer already writes `inputs.json` and `artifacts.json`. Bind
the summary and input inventory through those artifacts, verify actual input
hashes against the current builder arguments, and check 36 unique variant
records. Keep nonestimable results as records; do not require favorable
arithmetic agreement. A synthetic stale-input or missing-variant control should
fail before an output directory is created. This request does not change the
scientific endpoint or previous protocols.

## Checks that pass

- The V3 and recovery protocol hashes, raw Table 11 workbook hash/size and schema
  binding are checked. The source export and type/missing/formula ledger are
  jointly bound through the export receipt. The archived source passes all
  13-member CRC validation; extracted design/activity member hashes match the
  old design and retained raw activity workbook.
- The frozen DNA panel and outcome-independent grouping are hash-checked. The
  builder reads only element IDs from the old panel; it never uses old ESS
  outcome fields to construct activity observations. All **55,276** required
  IDs exist once in the source aggregate. No source alias is guessed, duplicate
  IDs are rejected even outside the cohort, and missing required IDs block the
  entire build.
- The published aggregate has 56,981 unique IDs. Its six shared-string missing
  cells occur outside the fixed cohort. All 165,828 required cells are ordinary
  finite numeric source cells, with no formula or exceptional cell type. The
  ledger rejection prevents a numeric-looking string or cached formula from
  silently becoming an observed measurement.
- The builder selects only the frozen `all_cell_types_summary` activity headers.
  It performs a decimal-lexeme-to-binary-float conversion and deterministic ID
  sorting. It does not center, scale, project, reconstruct counts, invert ESS,
  impute, clip, threshold or select by outcome. Zero and negative activities are
  retained by the synthetic controls. The exact decimal lexemes remain in the
  source export for precision checks.
- The official archive and selected member identities are retained. Final
  Table 10 is byte-identical to the previously reference-validated design
  workbook; final Table 12 is byte-identical to old S9. These are source copies,
  not independent validation. The builder's activity observations come from
  Table 11, not from Table 12 or model predictions.

## Reuse and interpretation

The final article grants CC-BY-4.0, including supplementary material unless a
specific third-party credit line says otherwise. The official source workbook's
document properties and all eight text header sets contain no additional reuse
terms; no custom license/copyright property member was found. The inspected
metadata and its hashes are retained in `activity-binding-source-checks.json`.
This check did not search arbitrary outcome-cell annotations. Source attribution
and original bytes remain retained.

Literal published-aggregate binding does not certify its reconstruction from
raw counts, reproduce all producer normalization history or establish independent
biological replication. The cohort inherits the earlier panel's ascertainment
and availability restrictions, and the source library was selected using earlier
activity measurements. Those limits persist even though this builder performs no
new outcome-based selection. Replicate arithmetic and release relationships
belong to the separately frozen recovery audit; predictive fitting remains a
separate gated action.
