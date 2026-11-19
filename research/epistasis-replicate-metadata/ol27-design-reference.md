# OL27 design accession and sequence-reference lead

The official ENCODE GEO submission records `ENCFF715XTT` as an upstream input
of the OL27 barcode map `ENCFF823YWM`. The same record reports the barcode map's
MD5, which agrees with the existing inventory. The IGVF catalog identifies
`IGVFFI9985FTZH` as its reference for the ENCODE library `ENCFF715XTT` and
provides a download URL. This establishes the accession chain down to the
catalog reference. The sequence match for the 12 candidate collision windows
still requires the file itself.
[ENCODE submission record](https://github.com/ENCODE-DCC/geo-submission/blob/fc2d3daabfcf90655147ba2f9d11680f3a8d1045/experiments/ENCSR517VAN_modified.json) ·
[IGVF accession mapping](https://github.com/IGVF-DACC/igvf-catalog/blob/1b9eca1d081920a8d2a2814a4b13e5971a65b475/data/data_sources/data_sources_file_fileset.yaml#L784) ·
[IGVF reference configuration](https://github.com/IGVF-DACC/igvf-catalog/blob/1b9eca1d081920a8d2a2814a4b13e5971a65b475/data/data_sources/data_sources.yaml#L1819)

The IGVF adapter expects a tab-separated design file with `name` and
`sequence` columns and reads effect records separately. A separately located
converter identifies the ENCODE input as the OL27 design FASTA. Its design
output preserves bases in uppercase and emits the individual aliases from
compound headers. Its selection requires annotation, liftover and membership
in upstream processed tables. The catalog association does not identify that
converter as the producer of the specific IGVF object. A TSV match can support
sequence identity. Missing aliases in a selected derivative leave original
library membership unresolved.
[IGVF adapter](https://github.com/IGVF-DACC/igvf-catalog/blob/1b9eca1d081920a8d2a2814a4b13e5971a65b475/data/adapters/mpra_adapter.py) ·
[Converter source](https://github.com/rnsherpa/Tewhey_MPRA_to_IGVF_spec/blob/836308447953f6f81711753e39d2104de4b28869/scripts/process_mpra.py#L637)

For the join, start with the original FASTA where available. Check the file
accession, release, content type and transferred checksum against source
metadata. Compare the complete design sequence for each exact candidate alias
against its S18 construct in the recorded orientation. Retain mismatches and
conflicting aliases separately. Require all four verified constructs per
window and keep the barcode-map status as a separate field. Compound aliases
can share one barcode count and therefore cannot provide independent
measurements. The input and reporter sample pairing and usable replicate
coverage remain separate checks.
[Existing barcode inventory](ol27-name-overlap.json)

Both publisher portals returned HTTP 200 and confirmed the files are released.
The original FASTA is a gzip file of 2,057,984 bytes with publisher MD5
`42c6e1a624a869e1f06656d84036a608`. The IGVF portal describes its public
2,888,359-byte file as MPRA sequence designs on GRCh38, with MD5
`052aca8a7c5cf0cad6521a2a6d39850f`. The retained metadata includes the
uncompressed MD5 for each and the URLs. The catalog API probe returned
HTTP 404 separately from these two successful publisher responses.
[Published metadata receipt](ol27-source-metadata.json) ·
[Metadata replay](https://github.com/mingley/helix/actions/runs/44611075449)

The crosswalk replay validates both checksums for each design file, recreates the
12 candidates from the pinned design and barcode map, and compares each of the
48 construct aliases to its audited S18 insert. It records distinct sequence
hashes under each alias, preserving missing, conflicting and mismatched hits.
It treats a semicolon-delimited member of a parenthesized FASTA header as
secondary alias evidence and reports the match type. The IGVF parser checks
the expected design header before comparing aliases.

```bash
python3 research/epistasis-replicate-metadata/audit_design_crosswalk.py --download
```

The downloaded inputs are temporary. A cached copy of the two design gzips
and the pinned OL27 barcode map can be replayed with `--inputs DIRECTORY`.

The source probe requests only object metadata and records response hashes,
published file hashes and provenance fields when available. Run it with
`python3 research/epistasis-replicate-metadata/inspect_design_reference.py`.
The GitHub workflow runs this probe when its source changes and also supports
manual replay. HTTP failures are reported per endpoint. The file and its
schema can be checked directly after obtaining the object. The project retains
the earlier outcome exposure history. The catalog adapter documentation also
displayed an example effect row for `HepG2_DNasePeakNoPromoter1` during this
review. Its values were unused.
