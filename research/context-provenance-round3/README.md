# Final published joint-library source discovery

The final Nature article links an accessible publisher archive containing
Supplementary Tables 1–13. Its exact URL is in [candidates.json](candidates.json)
and [publisher-attachment-head.json](publisher-attachment-head.json). The URL was
extracted from the public Nature article's attachment link; it was not constructed
to evade a challenge on another service. HTTP HEAD returned 200, a content length
of **342,234,089 bytes**, byte-range support and a last-modified timestamp of
14 February 2025. The ETag is an opaque validator, **not a verified MD5**.

This discovery round acquired no supplementary archive bytes, member data,
workbook headers or quantitative outcomes. Acquisition and endpoint inspection
require the centrally frozen round 3 protocol. No model was fit and no author
was contacted. Parent owns the primary acquisition and parser.

The final article distinguishes these tables:

| Final table | Article evidence | What remains unknown |
| --- | --- | --- |
| 10 | Joint-library composition/design, and the full element specificity score table (paragraphs Par26 and Par57) | Member filename, schema, and whether upstream activity columns are included |
| 12 | Same-element assignments across cell types to the same cross-validation fold (Par50) | Whether it contains upstream observed activity, ESS, predictions or only fold assignments |
| 8 | Large-scale library folds and model performance | Not established as a joint-library raw-activity source |
| 9 | Variant-effect/allele-specific validation contingency data | Different role from old Zenodo's Table S9 |

Final table numbering therefore differs from the older Zenodo deposition.
Published Table 10 should be identified by the final archive inventory, not by
assuming that older workbook names or labels are unchanged.

The Methods define the published ESS as a within-cell z-score of activity minus
the element's mean z-score across all three cells. The upstream activity is the
mean over three independently infected replicates after each replicate's
log2(RNA/DNA) values are centered on that replicate's full-library median. These
endpoints have different scales. An ESS table alone does not recover the cell
offsets, scales or original activity. The protocol must preserve the distinction
and fail closed if the required columns are absent or ambiguous. Our earlier raw
activity v1 semantic failure remains intact. Earlier rounds already exposed the
older Zenodo design and ESS panel; rediscovery of the same panel is not an
independent validation cohort.

The article is licensed CC-BY-4.0, including supplementary material unless a
specific third-party credit line states otherwise. The exact permissions block
is retained in [article-source-metadata.json](article-source-metadata.json) and
[final-article.xml](final-article.xml). Attachment-specific credit lines and terms
cannot be checked until acquisition. Cite the authors and DOI
[10.1038/s41586-024-08430-9](https://doi.org/10.1038/s41586-024-08430-9), preserve
the license and indicate transformations.

The PMC article also publishes the archive link and a rounded 326.4 MB size, but
its attachment HEAD returned an HTML challenge response. That is not evidence of
a downloaded ZIP. The older guessed static-content hostname failed CONNECT 403
and was abandoned. Normal official PMC OA metadata paths returned 404;
EuropePMC's supplementaryFiles HEAD timed out; its article-attachment route
returned 403. No access control was bypassed. Current runtime policy was
enforced unrestricted, and inherited proxy and TLS verification were preserved.

Receipts retain only relevant safe HTTP headers: representation, range,
cache/server/date and redirect metadata. Cookies, nonce-bearing fields and opaque
reporting headers are omitted. Full public Nature/PMC article HTML was used to
extract attachment metadata, hashed, then omitted; relevant link/caption extracts
and body hashes remain. The stable EuropePMC JATS article is retained byte for
byte; its SHA256 matches the earlier round's article hash. `manifest.json`
records the retained files and hashes. No shared frozen files were changed.
