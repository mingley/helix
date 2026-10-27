# Reconciled GRCh38 HBB reference

This bundle fixes one explicit reference model: human HBB on GRCh38 chromosome
11, RefSeq mRNA NM_000518.5 and protein NP_000509.1. It supports sequence-derived
mechanism calculations independently of the molecular-consequence labels in the
ClinVar demonstration fixture.

## Sources and integrity

The five input files are unchanged public HTTPS responses:

| File | Source | Purpose |
| --- | --- | --- |
| genomic.fasta | NCBI EFetch, NC_000011.10:5225001-5227500 | Forward-strand GRCh38 genomic bases |
| transcript.fasta | NCBI EFetch, NM_000518.5 | Versioned spliced mRNA |
| protein.fasta | NCBI EFetch, NP_000509.1 | Versioned 147-residue precursor translation |
| transcript.gb | NCBI EFetch, NM_000518.5 GenBank | CDS interval and protein accession |
| annotation.json | UCSC hg38 ncbiRefSeq track | Exon/CDS coordinates and strand |

The manifest records each exact request URL, retrieval time and SHA-256 digest.
UCSC exposes NCBI RefSeq annotation on the assembly; agreement between these
files is a cross-source representation check, not independent biological
replication. Accessions pin the sequence versions. Annotation metadata can change
even when the sequence accession version remains unchanged, so all raw responses
are also committed by their hashes.

Sources and terms:

- [NCBI RefSeq NM_000518.5](https://www.ncbi.nlm.nih.gov/nuccore/NM_000518.5)
- [NCBI RefSeq NP_000509.1](https://www.ncbi.nlm.nih.gov/protein/NP_000509.1)
- [UCSC API documentation](https://genome.ucsc.edu/goldenPath/help/api.html)
- [NCBI policies](https://www.ncbi.nlm.nih.gov/home/about/policies/)
- [UCSC download terms](https://genome.ucsc.edu/license/)

Source data retain their original terms; the repository code license does not
relicense NCBI/UCSC content.

## Coordinate model

All manifest coordinates are zero-based half-open intervals. The genomic
sequence covers [5225000, 5227500). Its FASTA header uses NCBI's one-based
inclusive notation.

HBB is on the **negative strand**. Exons in the manifest are ascending genomic
intervals, which is the reverse of their order in the mRNA:

| Genomic interval | Length | Position in transcript |
| --- | --- | --- |
| [5226929, 5227071) | 142 nt | First exon |
| [5226576, 5226799) | 223 nt | Second exon |
| [5225463, 5225726) | 263 nt | Third exon |

The transcript spans [5225463, 5227071); genomic CDS bounds are
[5225597, 5227021). These CDS bounds include intervening introns. The actual CDS
is obtained by intersecting each exon with those bounds, splicing, then reverse
complementing. Translating the entire genomic interval between CDS bounds would
be incorrect.

The first and second introns in transcript order are 130 and 850 nt. Both have
GT donor and AG acceptor dinucleotides in transcript orientation. A coordinate
near an exon boundary must be interpreted with strand and transcript context.
These motifs alone do not establish whether a variant disrupts splicing.

## Verified reconciliation

The importer refuses to write a bundle unless all of these checks pass:

1. Exact versioned transcript, gene, primary chromosome, strand and complete CDS.
2. Valid ascending, non-overlapping exon intervals inside the genomic reference.
3. Reverse-complemented spliced genomic exons exactly equal the 628 nt RefSeq mRNA.
4. Spliced genomic CDS equals mRNA interval [50, 494), which GenBank describes as
   CDS 51..494 in one-based inclusive notation.
5. Standard genetic-code translation of the 444 nt CDS equals the 147-residue
   protein followed by one terminal stop, with the expected ATG start.

The protein retains the initiator methionine. Historical mature beta-globin
numbering often omits it: the familiar sickle substitution called beta6 Glu>Val
is p.Glu7Val against this RefSeq protein. Sequence calculations must not silently
mix those numbering conventions.

These are computational reference checks. They do not prove transcript abundance,
tissue-specific regulation, effect sizes, clinical pathogenicity, or a new
biological finding. The model covers one transcript and one reference allele.

## Use and refresh

Run the offline integrity and reconciliation check:

    python3 scripts/fetch_hbb_reference.py datasets/hbb-reference --verify

Fetch another snapshot into a new directory:

    python3 scripts/fetch_hbb_reference.py datasets/hbb-reference-new

Existing directories are refused. Downloads use normal certificate verification,
inherited network/proxy settings, bounded retries for transient errors, and
NCBI-compatible request pacing. No credentials or private sequences are needed.
The original responses are preserved byte-for-byte, including trailing blank
lines in NCBI FASTA files.
