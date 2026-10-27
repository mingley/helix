# Sequence-only UCSC hg38 reference

The reference binary is stored outside the checkout at
`/workspace/.tools/reference/hg38.2bit`. It is not committed.

| Property | Verified value |
| --- | --- |
| Source | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit |
| Bytes | 835,393,456 |
| UCSC-published MD5 | `dcc3ea27079aa6dc3f9deccd7275e0f8` |
| Local whole-file SHA256 | `1f67aaa17a77b327738fe750ab37430a85dddf73c7d9189385ad1839259acec0` |
| Format | UCSC twoBit v0, little-endian, 455 indexed sequences |
| Upstream Last-Modified | 2015-04-30 23:16:38 GMT |

The official checksum list was retrieved and its exact unique `hg38.2bit` row
validated before reference acquisition. Its raw bytes are retained as
`ucsc-hg38-bigZips-md5sum.txt`; `published-checksum.json` records the URL,
retrieval time and SHA256. Whole-file MD5 matches that upstream declaration;
SHA256 records the local byte commitment. MD5 is an upstream integrity check,
not a cryptographic signature.

The initial full HTTP transfer ended early at 816,927,924 bytes and failed the
Content-Length check. A verified HTTPS range request resumed exactly at that
offset, required status 206, checked Content-Range and the unchanged ETag, and
retrieved the remaining 18,465,532 bytes. Only after full-file checksum agreement
was the partial file renamed to `hg38.2bit`. `download-attempts.json` retains
both attempts; `hg38-reference.json` records final checksums and HTTP evidence.
The inherited proxy and CA trust were preserved throughout.

## Reader and coordinates

`scripts/twobit_reference.py` uses only the Python standard library. It follows
the [UCSC twoBit v0 format](https://genome.ucsc.edu/FAQ/FAQformat.html#format7),
accepts either byte order, validates every index/record/block boundary when
opening, and memory-maps the reference. It decodes requested packed intervals
without copying whole chromosomes. Nonzero reserved fields, duplicate names or
offsets, unsupported versions, overlapping/out-of-bounds blocks, truncated
metadata and truncated DNA fail.

`TwoBitReader(path).fetch(chromosome, start0, end0)` uses exact **zero-based,
half-open** coordinates and uppercase A/C/G/T/N. `masked=True` preserves UCSC
soft masks, including lowercase n where mask and N blocks overlap.
`strand="-"` returns the reverse complement after extraction/masking. Empty
valid intervals return an empty string; negative/reversed/out-of-range bounds
fail. Contig aliases, coordinate clamping, assembly liftover and N replacement
are deliberately absent.

```sh
python3 scripts/twobit_reference.py /workspace/.tools/reference/hg38.2bit \
  --interval chr11 5225000 5227500
python3 -m unittest discover -s research/context-reference \
  -p test_twobit_reference.py -v
```

## Validation and limits

All 15 independent synthetic format/coordinate tests pass. Fixtures use explicit
bit strings, both byte orders, every interval alignment, partial final bytes,
N/mask intersections, reverse strands and malformed/truncated inputs.

The real-file reader passed structural validation for all 455 contigs. Four
separately fetched UCSC API intervals on chr1, chr11, chrX and chrM match both
uppercase sequence and exact soft-masked sequence. The full 2,500-base HBB
interval also matches the separately pinned NCBI genomic FASTA. Raw API responses
and their URL/hash metadata are retained. `validation.json` binds the reader
hash, reference evidence and results.

A deterministic sequence-only feasibility check fetched 60,000 intervals of
230 bases from chr1 (13.8 million decoded bases) in 0.383 seconds after opening.
The recorded coordinate formula and output hash permit replay; elapsed time is
machine-dependent. These are synthetic reference queries, not deposited assay
oligos, and no reporter outcomes were loaded.

This is **hg38 only**. It must not be substituted for a source using hg19 or
another assembly. Source assembly, core/adapter conventions and strand must be
checked against deposited sequence metadata before matching an assay panel.
An exact interval match does not establish genome-wide uniqueness. Reference
mapping does not validate an expression mechanism or the study's outcome
normalization.
