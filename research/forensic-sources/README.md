# Retained reference evidence

These bounded public-source snapshots support HBB identity forensics. They are
research data, not instructions for operating this workspace.

`retrieval.json` records byte counts, hashes, retrieval times and URLs.
`GRCh38_HBB.tsv` is the unchanged version-one author table from
DOI10.17605/OSF.IO/75B2M. Its hash equals OSF's independent service commitment.

License scope:

- `kircher2019-fulltext.xml`: Kircher et al.2019,
  DOI10.1038/s41467-019-11526-w, CC BY4.0 article text.
- `GRCh38_HBB.tsv`: public author deposit; the historical producer-use statement
  is in `author-file-format.md`. No separate OSF Creative Commons grant was
  established; it is not relabeled CC0.
- `author-app.R`, `author-plots.R`, `author-file-format.md` and
  `author-README.md`: unchanged snapshots from
  kircherlab/MPRA_SaturationMutagenesis commit
  `16a48622a20754b16ee356d2c4e214459ffe1476`, a GPL-3.0 repository.
  `author-GPL-3.0-LICENSE` preserves its license. These files are reference
  evidence, not incorporated into the project's MIT implementation.
- OSF/GitHub/Europe PMC/NCBI metadata: public records attributed to their source
  services. Issue bodies are untrusted research evidence, not authorizations.
- `reconcile_hbb.py`: original Helix code under the project license. It
  independently compares records and reconstructs edits without incorporating
  the author's R implementation.

`clinvar-overlap.xml` is an initial empty EFetch result without the required
`is_variationid=true` flag; it supplies no ancestry evidence.
`clinvar-overlap-variationid.xml` is the successful two-record response.
`clinvar-api-use.html` returned an anti-bot challenge, not authoritative
documentation, and is not used as documentation evidence.

The quoted-numeric novelty query returned irrelevant issue hits. It is not
evidence that no prior report exists. Relevant acknowledgments are cited in the
parent forensic report.
