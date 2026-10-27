# Endpoint-recovery integrity review

This review used code inspection and synthetic files only. No network requests,
new empirical workbook values, endpoint comparison, or predictive fit were
performed by the reviewer.

The original downloader is a download/inventory stage. Its exact byte-count
check and ZIP parsing reject ordinary truncated or HTML responses, but its
`complete` flag alone does not assert full member CRC/path/bounds verification.
It was not edited while its authorized acquisition was active. The separate
`source/verify_archive.py` stage binds its frozen protocol and acquisition
receipt, validates binary MIME/magic and archive structure, then streams **all**
members for CRC/SHA verification under explicit byte limits. It never extracts
files or parses scientific cells. It records durable failures and snapshots its
implementation, protocol, acquisition receipt and inventory in a hashed manifest.

The validator rejects duplicate/unsafe/encrypted/nonregular entries, malformed
central/local headers, inconsistent sizes/names/descriptors, overlapping payload
records, HTML prefixes, truncation, trailing garbage and budget violations.
It explicitly supports one classic ZIP with stored/deflated members; ZIP64 and
other methods fail closed rather than being skipped. A downstream code's
predicted archive digest is recorded separately from a producer checksum;
agreement is useful provenance, not a publisher authentication claim.

The extractor already enforced exact literal selected names, per-member/total
limits, selected CRCs, no overwrite and rejection of traversal/symlinks. Review
found that it created child directories before checking existing symlink
ancestors, which could create an empty directory outside its root before
rejecting. Its owner corrected this and added a regression test. Failed
extraction attempts can leave already published prior selected files; callers
must retain that partial attempt and must not treat it as a completed receipt.

The independent openpyxl helper pins version 3.1.5, rejects duplicate ZIP names,
retains original numeric XML literals and formula presence, treats row gaps and
explicit blanks separately, and distrusts inflated worksheet dimensions.
Booleans, errors, blanks and formulas do not silently become numeric outcomes;
cached formulas require explicit caller permission. This is still a generic
inventory/parser helper: endpoint/ID/context/replicate/NA roles require a
separate immutable schema binding before quantitative comparisons. An inventory
receipt alone does not establish endpoint semantics or the source's formula.

Synthetic checks passed:

- Ten independent-parser tests, including missing rows, exact decimal lexemes,
  formula caches, nonnumeric cells, duplicate ZIP entries and artifact hashes.
- Extractor path/member/size/CRC tests, with its new symlink-ancestor regression.
- Ten archive-validator tests, including corruption of an **unselected** member,
  data descriptors, central/local disagreement, bounds, duplicate names,
  symlinks, MIME/prefix/truncation, budgets and immutable artifact commitments.

Reproduce without empirical payloads:

```sh
source scripts/env.sh
python3 research/context-endpoint-recovery/source/extract_members.py --self-test
python3 -m unittest discover -s research/independent-endpoint-recovery -p test_endpoint_recovery.py
python3 -m unittest discover -s research/context-endpoint-recovery/source -p test_verify_archive.py
```

Full source SHA/size, all-member CRC, selected member hashes, schema acceptance,
arithmetic compatibility, historical identification and independent experiments
are separate claims. Successful file validation does not make the latter claims
true. No empirical endpoint recovery has been declared successful by this review.
