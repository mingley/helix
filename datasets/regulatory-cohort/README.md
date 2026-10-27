# Metadata-selected human promoter cohort

This frozen CC0 MaveDB cohort contains nine primary Kircher2019 promoter elements,
five validation-only shared-locus contexts/library repeats, and previously exposed
HBB reference data. See [the dataset and selection report](../../research/regulatory-cohort.md)
for the exact inclusion policy, accession roles, mapping, units, missingness,
permission, limitations and reproduction commands.

`metadata-manifest.json` and `metadata-selection.json` commit to the metadata-only
phase. `frozen-protocol.json` and `outcome-acquisition-gate.json` record the exact
protocol and roles before new scores were requested. `manifest.json` commits to
all consumed inputs, raw sources, mappings and source-row accounting.

All 18,468 original rows are retained in the score CSVs. The strict SNV view has
15,699 rows; 2,769 source equality descriptions remain explicit and are not
reinterpreted as neutral controls or deletions. Numerical quality filtering and
prediction evaluation belong to the frozen protocol, not to ingestion.

A held-out promoter from this one study is not an independent study. Reporter
coefficients retain log2 effect units and do not become clinical probabilities.
