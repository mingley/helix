# Frozen primary-chromosome RefSeq context atlas

This dataset maps **curated annotation coverage**, not the complete human gene
repertoire and not the fraction of the human genome with biological function.

## Sources and scope

The raw files are exact HTTPS responses from UCSC's `hg38` database:

- `ncbiRefSeqCurated.txt.gz` and its SQL schema: curated genePred models.
- `ncbiRefSeqLink.txt.gz` and its SQL schema: NCBI Gene IDs, names and biotypes.
- `chromInfo.txt.gz` and its SQL schema: reference sequence lengths.

The file timestamps recorded in `manifest.json` identify these particular UCSC
snapshots; they are not asserted to be the latest available annotations. The
RefSeq tables here have an August 2025 last-modified date. Source URLs, retrieval
time, byte lengths, and SHA-256 commitments are preserved in the manifest.

The RefSeqLink table omits the 13 mitochondrial models present in the curated
track. `mitochondrial-reference.xml` is the exact NCBI EFetch GenBank XML response
for `NC_012920.1`. Its CDS protein accessions and `GeneID` cross-references resolve
those 13 models without guessed IDs. Their `protein_coding` category follows the
CDS features. The curated track does not contain the mitochondrial RNA genes;
this is an explicit source limitation, not evidence that those genes do not exist.

Only `chr1`–`chr22`, `chrX`, `chrY`, and `chrM` are retained. Alt loci, fix patches,
random sequences, and unplaced/unlocalized sequences are excluded. Counts of all
excluded contigs and models remain in `genes.json`. Primary chromosome lengths
are checked against the assembly-specific GRCh38 lengths in `genome-types`.

Sources and data terms:

- [UCSC hg38 database downloads](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/)
- [UCSC RefSeq track description](https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=refSeqComposite)
- [UCSC download/data-use guidance](https://genome.ucsc.edu/FAQ/FAQdownloads.html)
- [NCBI mitochondrial RefSeq NC_012920.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1)
- [NCBI policies](https://www.ncbi.nlm.nih.gov/home/about/policies/)

Raw source data retain their source terms; the repository's MIT code license
covers the importer and loader, not a new license on those source records.

## Normalization

All intervals are zero-based and half-open. Curated models are joined to their
NCBI Gene metadata by exact versioned accession; mitochondrial protein models
use exact protein accessions in the GenBank source. Missing or ambiguous mappings
fail the import. Exons must be ordered, nonoverlapping and inside the transcript;
CDS portions are intersections of exons with annotated CDS boundaries.

A gene-locus record groups `(NCBI Gene ID, chromosome, strand)`. It contains a gene
envelope and merged unions of transcript spans, exons and CDS portions across
its curated models. It does not pick one preferred isoform. Coverage is computed
from merged unions over all loci, so overlapping genes and isoforms contribute a
base only once. Gene envelopes are not substituted for transcript union coverage.

`genes.json` is this normalized catalog; `summary.json` is the derived report.
`genome-atlas::Atlas::open` verifies all commitments, coordinates, interval nesting,
counts and source metadata, recomputes the report in Rust, and compares it with
`summary.json`. An additional offline importer verification rederives the catalog
from the compressed raw sources, independently of the stored normalized file.

## Frozen snapshot results

| Quantity | Value |
| --- | ---: |
| Primary reference length, including mitochondrial DNA | 3,088,286,401 bp |
| Retained curated transcript/gene models | 93,455 |
| Gene/chromosome/strand loci | 29,949 |
| Distinct NCBI Gene IDs | 29,922 |
| Protein-coding loci | 19,427 |
| Noncoding RNA loci (explicit RNA biotypes) | 9,334 |
| Pseudogene loci | 1,188 |
| Other/unclassified biotype loci | 0 |
| Total loci outside protein_coding | 10,522 |
| Transcript-span union coverage | 1,421,275,766 bp (46.0215%) |
| Exon union coverage | 93,083,895 bp (3.01410%) |
| CDS union coverage | 34,561,429 bp (1.11911%) |
| Excluded models on non-primary contigs | 11,564 on 438 contigs |

Other-biotype loci include 6,716 lncRNA, 1,915 miRNA, 1,184 transcribed pseudogene,
541 snoRNA and other categories. Full counts are reported without collapsing all
non-protein-coding categories into a single biological meaning. The 27-locus
difference between loci and distinct Gene IDs reflects the grouping rule; gene
identity and placement are different quantities, including X/Y placements.

These counts describe this curated source set. They should not be compared with
GENCODE or Ensembl totals as if every source used the same gene definitions,
release, inclusion policy or unit of counting.

## Reproduce and inspect

The ordinary loader and these checks are offline:

```bash
source scripts/env.sh
cargo test -p genome-atlas --locked
python3 scripts/fetch_genome_atlas.py --self-test
python3 scripts/fetch_genome_atlas.py --verify datasets/genome-atlas
```

To fetch a new snapshot, choose a directory that does not exist:

```bash
python3 scripts/fetch_genome_atlas.py datasets/genome-atlas-new
```

The importer never overwrites an existing snapshot. Source requests are bounded,
TLS verification remains enabled, and no credentials or private data are needed.
