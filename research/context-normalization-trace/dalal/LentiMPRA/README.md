# LentiMPRA processing outputs

The build writes four Parquet files from Supplementary Tables S3, S4, S8,
S10, S11, and S12.

| File | Rows | Contents |
|---|---:|---|
| `lentimpra-hepg2.parquet` | 164,307 | HepG2 large-scale library |
| `lentimpra-k562.parquet` | 243,780 | K562 large-scale library |
| `lentimpra-wtc11.parquet` | 75,542 | WTC11 large-scale library |
| `lentimpra-joint.parquet` | 60,000 | Joint three-cell library |
| `published-performance-large-scale.parquet` | 180 | Published S8 per-fold Pearson correlations |
| `published-performance-joint.parquet` | 150 | Published S12 per-fold Pearson correlations |

The three cell-specific files share one 23-column schema. The joint file has a
44-column schema with separate HepG2, K562, and WTC11 measurements and one
published 10-fold assignment.

Every design row is retained. Missing assay, QC, or modeling values remain
null. `sequence` is the 200 nt assayed insert after removing the two constant
15 nt synthesis adapters and is the actual model input; reverse constructs
already contain the assayed reverse-complemented insert.

Columns are ordered as identity and sequence, hg19 then hg38 location,
orientation, numerical measurements and targets, binary activity, and fold.
The two assemblies require explicit `hg19_chrom`/`hg38_chrom` namespaces
rather than one ambiguous generic coordinate. Starts and ends are zero-based,
half-open, BED-like source coordinates; exceptional source interval widths are
preserved. The 345 non-200 nt hg38 intervals comprise 153 promoters, 44
potential enhancers, and 148 Inoue et al. controls; all hg19 intervals are
200 nt. Assembly-specific strand is source genomic strand, while
`orientation` is construct orientation.

The two performance files are long-form transcriptions of the paper's
`Benchmarking results` sheets. They contain `dataset`, `cell_line`,
`test_fold`, `model`, and `pearson_r`; no model outputs are recomputed.

Run:

```bash
bash scripts/create_environment.sh
bash scripts/run_all.sh
```

The root [`manifest.json`](../manifest.json) pins every source and output,
ordered Arrow schemas, environment and processing identity, assay-specific
targets, coordinate semantics, and the cards-plus-Parquets Hub allowlist.
Normal processing only installs outputs after they match that manifest.
The row-level source reconciliation and recovered transformations are
documented in [PROCESSING_NOTES.md](PROCESSING_NOTES.md).
