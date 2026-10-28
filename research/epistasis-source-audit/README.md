# Design audit: Siraj et al. variant-pair oligonucleotides

This is a **design-only inventory**, not a reporter-effect result or evidence of
biological epistasis. The source is Supplementary Table 18 of Siraj et al.,
*Nature* (2026), [doi:10.1038/s41586-026-10121-6](https://doi.org/10.1038/s41586-026-10121-6),
[publisher file](https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41586-026-10121-6/MediaObjects/41586_2026_10121_MOESM9_ESM.xlsx).
The unchanged XLSX is retained here as `siraj2026-S18.xlsx`; its download was
successful on 2026-10-28 UTC, it is 1,617,092 bytes, and its SHA-256 is
`92929b4cef82797237c4c2be770feb744daa97255ac6dbfe3d53d5a34c996fbb`.
The article states [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/),
subject to exceptions for separately credited material. The separate Zenodo
record is CC BY-NC 4.0 and is not the source of this XLSX; its terms are not
interchanged. The article and authors are the source; the audit transforms
their table into the derived JSON. See the prior
[source and license review](../variant-transfer/source/README.md).

The publisher labels S18 “Oligo sequences for variant pairs.” The sheet has
the seven fields `Variant 1`, `Variant 2`, `Distance`, `Center Variant`,
`Window`, `Diplotype`, and `Oligo`. The authors selected 2,522 pairs within
150 bp using fine-mapping and CRE-overlap criteria, and designed up to six
200-bp windows per pair; thus this is a selected source cohort, not a sample
of all human variants. The methods report at least five transfection replicates
for each library, with no more than two per day and parallel replicates from
independently expanded cell batches. This design table does not report
replicate counts or response values by construct.

## Results and eligibility

All 35,616 data rows make 8,904 unique ordered variant-pair/center/window
groups, representing all 2,522 pair identifiers. Each group has exactly one
each of `RR`, `AR`, `RA`, and `AA`; no labels are missing, duplicated, or
extra. There are 33,184 individual 200-base sequences and 2,432 shorter
sequences; all source sequences use A/C/G/T only.

The audit marks a group eligible if both IDs are distinct same-chromosome
primary-chromosome A/C/G/T single-nucleotide substitutions, their coordinate
distance is at least three and agrees with all four source distance fields,
all four sequences are literal 200-base A/C/G/T strings, and the labels form
an exact square. `RR` to `AR` changes precisely variant 1's stated reference
base to its stated alternative, `RR` to `RA` does the same for variant 2, the
two edit positions have the stated distance, and `AA` is exactly the combined
edit. This strict forward-allele interpretation is checked directly on the
source sequence; it does not rely on assuming a strand from the center label.
It yields **7,542 eligible design windows across 2,156 pair identifiers**.
Among the other groups, 1,170 have a non-SNV or other nonmatching identifier
and at least one non-200-base sequence; 192 have two SNV identifiers less
than three bases apart. These reasons can overlap, so the reason counts are
not a partition. Source labels and row numbers for all groups and the precise
reasons are retained in `s18-design-audit.json`.

The 7,542 eligible groups have 6,553 distinct reference-sequence SHA-256
fingerprints; 1,775 groups are in a nonunique-fingerprint class, and 757 such
classes span multiple pair identifiers. Hence windows and even different pair
IDs cannot be assumed independent. This is only exact reference-sequence
reuse, not a complete dependence or homology screen. It is not an exclusion
rule and does not infer which observations will be available.

## Exposure, replay and limitations

The companion S20 outcome workbook was previously downloaded as part of the
S30 bundle. Before this audit its S30 control sheet had been interpreted, but
S20 data rows had not been analyzed. During this audit its first two header
rows only were interpreted: per-window four-diplotype `Log2FC` fields are
present. To inspect S29 source metadata, `openpyxl` loaded the bundle's shared
string table, which could instantiate strings from other sheets, although no
S20 data row or outcome value was interpreted. This history and prior source
study exposure preclude describing this as a new, completely unexposed study.
Subsequent upstream-header inspection also buffered 63 bytes after the header
without printing them; the [endpoint review](../epistasis/source/siraj-s20-semantics.md)
records that separate exposure and its scope.
Methods identify the activity scale as a DESeq2 RNA-versus-DNA log2 fold change
after library-specific summit-shift normalization. The authors' pinned public code
[`fig3_file_overview.R`](https://github.com/julirsch/finemapped_mpra/blob/1e38d1e3d4414dfa3a62b22110eba0136eda197d/code/fig3/fig3_file_overview.R)
references `code/tables/table20.R`, but that file is absent from the pinned
repository tree; fetching that exact raw path returned HTTP 404. The adjacent
`haplos_preprocess.R` reads separate `all_windows.txt.gz` and
`meta_windows.txt.gz` inputs and joins them by pair and cell type.

The separately licensed Zenodo record's
[`code.zip`](https://zenodo.org/api/records/15297965/files/code.zip/content)
(SHA-256 `0232a4ab60d0fec6eb3f1f3ed1ba0cbe08920fca821647beafc2c0b665b9b87e`)
contains `code/tables/table19.R` (SHA-256
`7c8cf5a9dfbee240ec7fbe495b1afb69f58910d42ad67b89726daca5a8820697`),
which appears to produce the S20 schema despite the numbering difference. It
groups by variant pair, centered variant, window and cell type, selects one
row with the maximum sum of square roots of the four mean plasmid counts,
breaks ties by the first row, renames the four upstream `*_Log2FC` fields and
pivots the centered variant and window into column names matching S20. This is
evidence for a shared selected source row for the four activities. An upstream
header also identifies a single library field on each such row. However, its
`stable19.txt` output has not been shown byte-identical to publisher S20;
the workbook linkage, upstream normalization and per-construct replicate support
still need verification before treating the published columns as a valid
interaction endpoint. This source code was inspected for provenance and was
not incorporated into the audit. The Zenodo record states CC BY-NC 4.0; it is
distinct from the publisher's S18 file. S20
availability, missingness, significance and response values were not checked
or used here. The S21 table is a selected, outcome-recoded non-additive set
and must not silently define an unbiased cohort. The authors have already
studied non-additive effects in this same experiment.

Run from the repository root with the configured research environment:

```sh
source scripts/research-env.sh
python research/epistasis-source-audit/audit_s18.py \
  research/epistasis-source-audit/siraj2026-S18.xlsx /tmp/s18-design-audit.json
cmp research/epistasis-source-audit/s18-design-audit.json /tmp/s18-design-audit.json
```

The source hash, worksheet name and headers are checked before parsing.
The audit is intentionally limited to insert sequences and source IDs; it
does not verify genomic assembly, source exclusion from prior studies,
regulatory mechanism, or a model's predictive performance. In particular,
design eligibility is not the number of measured or independent interactions.
