# Siraj OL27: replicate-source and lookup feasibility

This is an **outcome-free source and design investigation**, not a frozen
interaction cohort or evidence of a biological effect. It identifies a
potential source of replicate-level barcode counts and explores whether a
published four-construct design might be represented in one library. It does
not establish the published S20 endpoint, sample pairing, usable construct
coverage, normalization, or a verified original-design ID-to-sequence crosswalk.
The paper, previous Siraj single-variant outcomes, and selected pair-design
information have already been explored in this project.

## The K562 source

The publisher's S29 data-source sheet identifies K562 library OL27 experiment
`ENCSR857JFG`, GEO series
[GSE266781](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266781),
and processed variant table `ENCFF190GIK`. A GEO sample README describes a
larger derived-from-FASTQ element-quantification table `ENCFF970OLE`, the
barcode lookup `ENCFF823YWM` among its inputs, and the source's uncompressed
size and MD5. A bounded gzip-prefix read returned exactly this first line:

```text
Barcode  Oligo  plasmid_r1  plasmid_r2  plasmid_r3  plasmid_r4  plasmid_r5  k562_r1  k562_r2  k562_r3  k562_r4  k562_r5
```

The actual separators are tabs. The compressed source is 212,873,667 bytes;
its first 65,536 bytes have SHA-256 recorded in [evidence.json](evidence.json).
The GEO README's uncompressed-size and MD5 claims are recorded but **have not
been checked against the full count table**. No numeric RNA or reporter outcome
row was decoded or interpreted. The range contained compressed bytes that may
encode such rows; the decompressor returned no byte past the header newline.

GEO series metadata independently list five K562 samples
`GSM8254484`–`GSM8254488`, with biological-replicate labels 1–5 and technical
replicate 1. The separate OL27 plasmid sequencing series
[GSE266843](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266843)
lists five cloning samples with biological-replicate labels 1–5. These labels
do **not** prove that count-column suffixes map to particular GEO samples, that
RNA and DNA indices form matched transfections, or that the cloning material
constitutes independent biological cultures. The article describes independent
transfection cultures and a constructed plasmid library; shared preparation and
systematic construct/barcode effects may persist across samples.

Published [MPRASuite source](https://github.com/tewhey-lab/MPRASuite/blob/e853506aaa1edbae395c5d0870c141754a8fe667/MPRAmodel/MPRAmodel.R)
aggregates barcode counts by oligo and uses condition information and DESeq2
summit-shift normalization. This historical version is consistent with the
source schema; it is **not authenticated as the exact code or configuration**
that produced this GEO file. Raw counts do not automatically reproduce S20,
and a contrast using a common estimated DNA denominator can induce shared
cross-replicate error. Replicate indexing, filtering, pairing, normalization,
sequencing depth, shared library effects and a variance model all require
source-specific resolution before confirmatory inference.

## Lookup-name inventory

S29 also names OL27 barcode lookup `ENCFF823YWM`. Its GEO source is available
under [sample GSM8254614](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8254614).
The downloaded 296,690,753-byte gzip has SHA-256
`58f0c48c6d0a21b4a558794f31dcf38ccd6e8c766d56450512a28a27de059a15`;
the complete decompressed file is 1,039,687,724 bytes and matches the GEO
README's MD5 `c18ba6775c64fbef1da73cede0d29dea`. The compressed source is
cached outside Git. The GEO metadata list `ENCFF715XTT` among its upstream
inputs without identifying that input's type; we could not establish an
accessible original oligo-ID-to-sequence manifest from that metadata.

The ten-field lookup has no header. The current
[MPRAmatch output documentation](https://github.com/tewhey-lab/MPRASuite/blob/da5861938be0b7f3f98539a4e029a3649b5a58f2/MPRAmatch/output_file_explanations.md)
describes field 2 as comma-separated oligo names and field 5 flags 0, 1 and 2
as single, conflicting and failed mappings, respectively. It has not been
proven to be the exact schema version used for this file. Our exploratory code
generates candidate strings from each S18 pair, reference/alternate labels,
center and window, using the visible naming pattern; it does **not** establish
that the strings represent the S18 DNA sequences. Among 7,542 design windows:

| Candidate-name count present in lookup | Windows |
| --- | ---: |
| Four | 256 |
| Two | 252 |
| Zero | 7,034 |

The 256 represent 105 pair identifiers. All four names for each of those
windows also appeared on at least one single-name, flag-0 row. Of the 282
previously identified feature-collision windows, **12 windows in five pair
identifiers** have all four candidate strings under that flag rule. These are
**syntactic overlaps only**. They cannot be used as a measured-outcome cohort
or as evidence that a lookup token matches a particular S18 sequence.
The [independent reconstruction](../epistasis-replicate-independent/README.md)
reproduces the inventory counts. Different windows can share constructs and
are not independent experiments.

The producer's `code/tables/table17.R` in the separately licensed Zenodo code
archive reads an internal `haplotypes_info.txt.gz` containing the design `ID`
and oligo fields, but omits the ID from its S18 export. That internal design
file was not present in the inspected Zenodo `data_preprocess.zip`, and the
inspected `resources.zip` contains a chain file only. Filename searches do not
prove that no source elsewhere can supply this crosswalk.

## Partial long-ID and sequence check

The separately licensed CC BY-NC 4.0 [Zenodo record 15297965](https://zenodo.org/records/15297965)'s
[`data_satmut.zip`](https://zenodo.org/api/records/15297965/files/data_satmut.zip/content)
is 93,320,206 bytes and contains the design member
`data/satmut/satmut_sequences.txt.gz`, 3,455,761 bytes with SHA-256
`dd0080c8d451710e766e6b2dadb35e952da680183d47eb5775adc7fbaf91dea0`.
It was obtained by bounded ZIP range reads; the member's ZIP CRC and its
complete SHA-256 were checked, but the full archive was not acquired or
independently hashed. Zenodo lists the full archive MD5 as
`0043848d7658ac5cf402634b9685102e`; the selected member occupies 2,507,667
compressed ZIP bytes and has ZIP CRC-32 `2316519572`. Satmut targets were selected in part using
previous assay outcomes; this is an **exposed, selected** design source. No
satmut response member was opened.

The [crosscheck ledger](satmut-crosscheck/pair-backgrounds.jsonl) retains all
148 `:m0` pair-background rows, including unmatching or ambiguous IDs and
non-200-base sequences. Candidate long identifiers are formed from S18's two
full variant identities, centers, windows and allele labels. The accompanying
[summary](satmut-crosscheck/summary.json) records exact sequence comparisons,
unmatched and ambiguous counts, and their relation to the OL27 syntactic names
and feature classifications. **124 of 148** pair-background rows match one
generated S18 long ID, covering 33 windows; all 124 sequences agree exactly.
The remaining 24 rows match no S18 long ID, none match multiple IDs, and six of
the 148 are not literal 200-base A/C/G/T strings. Four of the 124 rows, in two
non-collision windows, have an inferred short name present on a flag-0 single
name row. None of the 12 syntactic collision windows is validated by this
selected source. Long-ID/sequence agreement does **not** prove that removing
the second variant's alleles and changing separators produces the correct short
identifier in the original OL27 lookup.

## Exposure and reproduction

For S29 metadata, the workbook reader loaded the workbook-wide shared-string
table; we interpreted S29 source metadata, not S20 data rows. The first-line
probe of the headerless lookup decoded and printed **one full mapping record**,
including barcode, oligo identifiers and QC/alignment fields. That preview is
not retained here. The subsequent full lookup inventory instantiated every raw
row including barcode and QC/alignment bytes, and interpreted only the oligo
names and overall mapping flag. No RNA or reporter-response table was read.
The count-table prefix was buffered as compressed bytes and only its header was
returned by the decompressor. The independent check reads the same lookup, not
an independent biological experiment.

The retained GEO series/sample metadata and their URLs and hashes are in
[`sources/manifest.json`](sources/manifest.json). From the repository root:

```sh
source scripts/research-env.sh
mkdir -p "$HELIX_TOOLS_DIR/epistasis-replicate-metadata"
# Keep any existing lookup. inventory.py will reject a hash mismatch.
if [ ! -e "$HELIX_TOOLS_DIR/epistasis-replicate-metadata/ENCFF823YWM.tsv.gz" ]; then
  curl --fail --location \
    'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8254nnn/GSM8254614/suppl/GSM8254614_ENCFF823YWM_element_barcode_mapping_.tsv.gz' \
    --output "$HELIX_TOOLS_DIR/epistasis-replicate-metadata/ENCFF823YWM.tsv.gz.download" &&
  printf '%s  %s\n' \
    '58f0c48c6d0a21b4a558794f31dcf38ccd6e8c766d56450512a28a27de059a15' \
    "$HELIX_TOOLS_DIR/epistasis-replicate-metadata/ENCFF823YWM.tsv.gz.download" | sha256sum --check - &&
  mv "$HELIX_TOOLS_DIR/epistasis-replicate-metadata/ENCFF823YWM.tsv.gz.download" \
    "$HELIX_TOOLS_DIR/epistasis-replicate-metadata/ENCFF823YWM.tsv.gz"
fi
# Fetch only the design member through bounded ZIP ranges, checking its ZIP
# metadata, CRC, complete member size and SHA-256. Does not overwrite a file.
if [ ! -e "$HELIX_TOOLS_DIR/epistasis-replicate-metadata/satmut_sequences.txt.gz" ]; then
  python research/epistasis-replicate-metadata/fetch_satmut_member.py \
    --out "$HELIX_TOOLS_DIR/epistasis-replicate-metadata/satmut_sequences.txt.gz"
fi
python -m unittest discover -s research/epistasis-replicate-metadata -p 'test_*.py'
python research/epistasis-replicate-metadata/verify_metadata.py
# Optional: publisher S29 row and a fresh bounded count-header request.
python research/epistasis-replicate-metadata/verify_metadata.py \
  --workbook "$HELIX_TOOLS_DIR/variant-transfer-v1/S30-bundle.xlsx" --live-header
python research/epistasis-replicate-metadata/inventory.py \
  --lookup "$HELIX_TOOLS_DIR/epistasis-replicate-metadata/ENCFF823YWM.tsv.gz" \
  --out /tmp/ol27-name-overlap.json
cmp research/epistasis-replicate-metadata/ol27-name-overlap.json /tmp/ol27-name-overlap.json
python research/epistasis-replicate-metadata/satmut_crosscheck.py \
  --satmut "$HELIX_TOOLS_DIR/epistasis-replicate-metadata/satmut_sequences.txt.gz" \
  --lookup "$HELIX_TOOLS_DIR/epistasis-replicate-metadata/ENCFF823YWM.tsv.gz" \
  --out /tmp/ol27-satmut-crosscheck
cmp research/epistasis-replicate-metadata/satmut-crosscheck/pair-backgrounds.jsonl /tmp/ol27-satmut-crosscheck/pair-backgrounds.jsonl
cmp research/epistasis-replicate-metadata/satmut-crosscheck/summary.json /tmp/ol27-satmut-crosscheck/summary.json
```

The lookup and satmut member are external to Git; their URLs, sizes and hashes
are recorded in [evidence.json](evidence.json) and above. The inventory checks
the full lookup SHA-256, uncompressed MD5 and size before writing a new output;
the satmut check verifies the lookup and member hashes. Output creation is
exclusive. The range reader checks every server `Content-Range`, the 93,320,206
byte archive length from the source, ZIP member sizes and CRC, and complete
member SHA-256. It does not download the entire ZIP, so it does not verify the
archive's published MD5. Downloading and parsing the complete count/response
table remains a separate, currently unperformed step.
