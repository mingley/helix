# Context literature evidence

This directory supports [the context interpretation review](../context-literature.md).
It contains minimal bibliographic metadata and original paraphrased notes;
it does not contain row-level assay outcomes, sequencing counts, supplementary
outcome workbooks, model arrays or redistributed article bodies.

- `publications.json`: eight citation records, including two corrections,
  projected from the official Europe PMC REST API. Abstracts are omitted.
- `review.json`: thirteen source-backed design/prior-art notes and four proposed,
  unexecuted questions, with interpretation limits and missing evidence.
- `retrieval-manifest.json`: source URLs, article-reading response hashes,
  reuse notes, failed/nonarticle responses, and hashes for the retained JSON
  projections. Article text was read in temporary storage and is not copied here.

The reading receipts distinguish access from permission. Agarwal 2025 and
Kircher 2019 state CC BY4.0; Inoue 2017 states CC BY-NC4.0. Public access to
Klein, Tewhey and Ulirsch through the official NCBI BioC service is not treated
as a broad redistribution grant. The repository's software license does not
relicense those articles or any candidate source data.

Official citation API:
<https://www.ebi.ac.uk/europepmc/webservices/rest/search>

Query: `EXT_ID:(27831498 OR 29717003 OR 33046894 OR 27259154 OR 27259153 OR 29474912 OR 39814889 OR 31395865) AND SRC:MED`,
with `format=json`, `resultType=core`, and `pageSize=50`.
The provider's license labels remain metadata; they are not independent grants.

Published methods and findings were reviewed. No newly acquired quantitative
activity data were inspected for this review. Outcome-informed selection of
the source joint library and full-library endpoint normalization are disclosed
in the narrative. This is not a claim that the research history was blind.

Two corrections could be identified bibliographically but their text/scope
could not be verified through the available official responses. No inference
about wholesale invalidity or harmlessness is made. Resolve correction scope
before quantitatively reusing the corresponding alternative study.

The final Agarwal paper calls its joint library/ESS table Supplementary Table10.
The Zenodo filename `Table S7 - joint lib design.xlsx` uses different numbering;
its title alone does not certify metadata-only content or identify it with the
final table. No workbook was fetched by this review.
