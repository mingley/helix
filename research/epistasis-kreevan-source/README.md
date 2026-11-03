# Candidate independent source: Kreevan et al. variant-pair MPRA

This is a **source-feasibility review**, not a frozen study protocol, an
eligible-sequence inventory, or a biological result. Kreevan et al.,
[*Massively parallel characterization reveals context-dependent and
non-additive regulatory effects of closely spaced variant pairs*](https://doi.org/10.64898/2026.09.10.750584),
posted a preprint on 11 September 2026. It reports a separate human reporter
experiment with all four ancestral/derived haplotypes for 7,285 SNV pairs. The
research question and aggregate findings in its abstract are already published;
this source is prior art, not a blank-slate claim of discovering epistasis. It
is separate from the previously studied Agarwal and Siraj experiments, but
sequence or locus overlap with those studies has not been assessed.

The [GEO series GSE349490](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE349490)
was public on 9 October 2026. Its sample descriptions and the manuscript's
methods report paired RNA and plasmid-DNA preparations for seven replicates in
each of K562 and GM12878. GEO describes a barcode count table with separate
columns for every replicate and material. The manuscript describes splitting
cells from a transfection for RNA and DNA extraction. The metadata do not by
themselves establish independently expanded cultures, independence of errors,
or absence of shared plasmid-library, barcode, or batch biases. These are
transfection/cell-culture measurements, not donor-level replication.

## Designs, source scale and current blocker

The methods specify **170 bp of genomic insert**, flanked by 15 bp adapters at
each end, and pairs separated by at most 75 bp. The sequence-selection metadata
include eQTL and enhancer filters; this is a selected cohort, not a random set of
human variants. The four labels use ancestral/derived alleles, which must not
be silently equated with reference/alternate alleles. GEO's sample processing
metadata specify hg19. The existing model's theoretical locality property
extends to fixed-length 170-bp inserts, but its fitted activity calibration was
trained on 200-bp inserts and does not automatically transfer.

The paper cites **Table S1** for the design sequences, but the reviewed
[bioRxiv supplement page](https://www.biorxiv.org/content/10.64898/2026.09.10.750584v1.supplementary-material)
lists one supplementary-figures PDF. The reviewed
[GEO supplement directory](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE349nnn/GSE349490/suppl/)
lists an S2 barcode-count file and S4/S8 cell-line results, with no separately
listed S1 design file. The virtual sequencing platform GPL34284 reports zero
design rows. This is a bounded search, not proof that S1 does not exist
elsewhere.

The S2 header has six metadata columns and 28 RNA/DNA sample columns. A
separate, tightly bounded inspection of the first data record
found that its `Oligo` field is `wide85_6702_ref1_alt2_T/A_anc/der`, a label,
not a literal inserted sequence. This one record does not establish the format
of every record. The naming convention, genomic interval and actual four
insert sequences remain unverified. A list of such identifiers would not by
itself establish exact SNV identities, strand, reference sequence or a valid
four-construct square. Scanning the complete compressed count file to inventory
those labels would require processing through response/count rows and is not
part of this review. No counts or response values were decoded or analyzed.

## Evidence and exposure

The [header capture](header.json), [first-field capture](first-oligo.json) and
[metadata audit](metadata.json) preserve observed columns, 28 GEO sample/SRA
identities, input SHA-256 values, and explicit unmet gates. The count file was
not downloaded in full and its full-file hash is unknown. For each bounded
capture, an 8,192-byte compressed HTTP range was requested; socket layers may
have buffered that range. For the header, 266 compressed bytes were fed to the
decoder and **zero** decompressed bytes beyond the newline were instantiated.
For the first record, 318 compressed bytes were fed; decompression was capped at
one output byte and stopped at the second tab. It instantiated the leading
15-base barcode and the `Oligo` field, with **zero** bytes decoded after that
tab. No alignment or count field was decoded. Prefix hashes in the two captures
do not authenticate the entire S2 file.

Search results displayed the paper's **aggregate abstract outcomes** before
this review. The article XML, including other sections, was downloaded and
loaded into memory, while the reproducible extractor uses only the specified
methods sections and license metadata. GEO's published processing text,
including aggregate coverage/QC descriptions, was inspected as metadata. No
individual-pair outcome table, response value, or count value was interpreted;
this is not a wholly outcome-unexposed study. S4 and S8 were not downloaded.

The preprint's XML identifies **CC BY 4.0** for the article; this review does
not assume that article terms automatically settle every separate GEO/SRA
file's reuse conditions. Exact external source files are hash-bound in the
audit. The article and sample metadata are not stored here; re-acquisition
requires the same bytes, and the audit stops if the public copies have changed.

## Reproduction

The following fetches the metadata and reads only the bounded S2 prefix described
above. The article XML contains published results even though the extractor
selects methods only. `capture-*` and `inspect` require new output files.

```bash
kreevan_dir=$(mktemp -d)
curl --fail --silent --show-error --location \
  'https://www.biorxiv.org/content/early/2026/09/11/2026.09.10.750584.source.xml' \
  --output "$kreevan_dir/article.xml"
curl --fail --silent --show-error --location \
  'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE349490&targ=self&form=text&view=full' \
  --output "$kreevan_dir/series.txt"
curl --fail --silent --show-error --location \
  'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE349490&targ=gsm&form=text&view=full' \
  --output "$kreevan_dir/samples.txt"
python3 research/epistasis-kreevan-source/inspect_metadata.py capture-header --out "$kreevan_dir/header.json"
python3 research/epistasis-kreevan-source/inspect_metadata.py capture-first-oligo --out "$kreevan_dir/first-oligo.json"
cmp research/epistasis-kreevan-source/header.json "$kreevan_dir/header.json"
cmp research/epistasis-kreevan-source/first-oligo.json "$kreevan_dir/first-oligo.json"
python3 research/epistasis-kreevan-source/inspect_metadata.py inspect \
  --series "$kreevan_dir/series.txt" --samples "$kreevan_dir/samples.txt" \
  --article "$kreevan_dir/article.xml" --header "$kreevan_dir/header.json" \
  --first-oligo "$kreevan_dir/first-oligo.json" --out "$kreevan_dir/metadata.json"
cmp research/epistasis-kreevan-source/metadata.json "$kreevan_dir/metadata.json"
python3 -m unittest discover -s research/epistasis-kreevan-source -p 'test_*.py'
```

The script rejects hash, sample-pairing, assembly, methods and header
discrepancies with explicit exceptions, including under `python -O`. Five
synthetic controls check complete sample pairing, wrong sample/assembly
identities, bounded first-field decoding, and malformed/truncated input. The
sample metadata and paper are the source of those checks; a code replay is not
independent experimental replication.

## Gates before any interaction outcome analysis

1. Locate a source-bound design-to-sequence mapping; validate all four complete
   170-bp inserts, the two exact substitutions and their spacing/orientation,
   genomic reference/ancestral allele distinctions, and source IDs. Do not infer
   a sequence from a label.
2. Before accessing count/response rows, freeze grouping, locus/sequence
   dependence and prior-exposure rules; source sample identities; low-support
   and missingness handling; normalization and pseudocount choices; contrast
   orientation; weights; uncertainty; and stop conditions. The article's
   source QC and high-confidence selection use observed data and must not be
   silently repurposed as an outcome-blind cohort.
3. Establish whether the paired RNA and plasmid-DNA extractions support the
   intended replicate estimand and account for common barcode/library and batch
   effects. Independent sequencing is not independent biological replication.
4. Interpret any finding on the specified reporter and activity scale; seek
   further discriminating evidence before claiming an endogenous mechanism.
