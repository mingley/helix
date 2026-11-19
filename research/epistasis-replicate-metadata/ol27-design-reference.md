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

The source probe requests only object metadata and records response hashes,
published file hashes and provenance fields when available. Run it with
`python3 research/epistasis-replicate-metadata/inspect_design_reference.py`.
The GitHub workflow runs this probe when its source changes and also supports
manual replay. HTTP failures are reported per endpoint. The file and its
schema can be checked directly after obtaining the object. The project retains
the earlier outcome exposure history. The catalog adapter documentation also
displayed an example effect row for `HepG2_DNasePeakNoPromoter1` during this
review. Its values were unused.
