# Prior-exposure DNA and locus inventory

All **55,276 Agarwal cores**, **15 Kircher assays at ten promoter loci**, and the
explored HBB genomic/ClinVar reference intervals are retained for external-study
exclusion. The ten Kircher loci comprise the nine previous non-HBB primary
promoters plus HBB; replicate, time and cell-context assays remain linked to
the same loci. This is an inventory of explored evidence, not a new cohort.

The complete [inventory manifest](inventory-v1/manifest.json) has SHA-256
`861d838582ce1f0d52d846e6d20ee5d6c0f28c059db257032c9cff7a1fae6faf`.
Its twelve listed files total 23,176,069 bytes. All thirteen files including the
manifest reproduce byte-for-byte under a new output directory. All seven
synthetic/inventory controls pass; [verification.json](verification.json) retains
that check and the pre-output manifest-format correction.

## Reference and molecular identity checks

Every Agarwal core and every full Kircher target matches the locked UCSC hg38
reference exactly at its declared interval and strand. The Agarwal union of
five training sets contains exactly all 55,276 retained panel IDs. Kircher's
source coordinates explicitly declare hg38 and reconcile as one-based inclusive;
this inventory converts them to **zero-based half-open GRCh38 intervals**.
No hg19 coordinates are relabeled as hg38, and no chain conversion is needed
for these old datasets. Exact matching within a reported interval does not prove
genome-wide uniqueness.

All **18,468 previously explored Kircher MaveDB rows** have a bijective
chromosome/position/reference/alternate match to the original author's
checksum-verified OSF tables: **15,699 substitution rows and 2,769 deletion
rows**, without applying a support or outcome filter. Repeated contexts produce
**11,007 distinct substitution edits and 1,744 distinct deletion edits**.
`n.position=` is not interpreted as neutrality: the original table's explicit
`Alt=-` establishes a deletion for each such row. This is a molecular-key
bijection; numeric scores, counts and historical export mechanics are not
verified by this inventory.

Two source details are retained explicitly:

- The author filename `GRCh38_GP1BA.tsv` contains chr22 alleles matching the
  GP1BB target. Full target/reference and every allele establish the mapping;
  the filename alone does not.
- The original MSMB table contains seven rows whose genomic identity fields are
  all `NA`. They are absent from the previously explored MaveDB table and are
  recorded separately. No DNA identity is invented for them.

The original assay uses additive coefficient estimates from mutagenized
libraries. Names such as `native-edits` mean isolated single edits applied to
the full source-length target, as used by the prior benchmark. They do not
establish that each isolated edited sequence was measured as a single-mutant
barcode construct, or resolve its original multi-mutant background.

## Files and exact homology representation

| File | Content |
| --- | --- |
| `agarwal-cores.jsonl.gz` | All 55,276 source-selected 200 nt cores, coordinates, strands, hashes and original training-fold membership |
| `kircher-loci.json` | All ten complete promoter references and all fifteen assay roles |
| `kircher-native-edits.jsonl.gz` | Every distinct source-resolved isolated substitution/deletion sequence, parent locus and all supporting identity records |
| `kircher-reference-windows200.jsonl.gz` | 1,694 derived reference windows with the observed single-substitution alternatives factorized by local position |
| `kircher-deletion-windows200.jsonl.gz` | 128,641 distinct canonical deletion-created windows absent from reference windows, with every edit/window occurrence |
| `hbb-auxiliary-loci.json` | The explored HBB gene-reference and ClinVar-reference intervals and DNA |
| `hbb-auxiliary-windows200.jsonl.gz` | All 3,601 derived 200 nt windows from those two auxiliary reference intervals |
| `molecular-identity-bindings.json` | Per-assay old-source identity bijections and unmatched author rows |
| `inputs.json` | Hashes of every retained old panel/cohort source plus exact reference and author receipts |
| `summary.json` | Coverage, reference checks, conventions and limits |
| `implementation.py`, `twobit_reference.py` | Exact builder and reference-reader snapshots |

Every DNA record includes its literal sequence, SHA-256 and canonical
reverse-complement SHA-256. Derived windows are exclusion substrings, not newly
asserted physical reporter constructs. HBB's 187 nt promoter and 186 nt deletion
sequences have no 200 nt windows; their actual genomic locus remains covered by
the declared coordinate guard.

The substitution factorization is exact. If a query differs from a reference
200 nt window at at most twenty positions, it already qualifies at 90% identity.
At twenty-one mismatches, an exposed substitution qualifies if and only if its
documented alternate matches the query at one differing position. More than
twenty-one mismatches cannot be rescued by one substitution. Use both strands,
reversing/complementing mutation positions and alleles with the window. An index
for this factorization must retain candidates up to twenty-one mismatches;
twenty-one disjoint exact segments alone only guarantee recovery up to twenty
mismatches. Use twenty-two segments, or another exhaustive method, for the
factorized branch. Deletions require the separately retained explicit windows.

## Reproduction and attribution

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover \
  -s research/variant-transfer/prior-exposure-inventory -p test_inventory.py -v
python3 research/variant-transfer/prior-exposure-inventory/build_inventory.py \
  --out /tmp/prior-exposure-review \
  --author-cache research/variant-transfer/prior-exposure-inventory/raw-author-tables
```

Choose a fresh output path. The builder uses the already locked hg38 reference
and retained old author tables, hashes every old source before interpretation and
never opens a Siraj file. The acquisition helper can recreate the old author
cache using only metadata-committed URLs/SHA-256s, normal TLS and bounded requests;
use `--receipt /tmp/old-author-acquisition-review.json` to preserve the original
receipt, and `--cache PATH` when selecting a different cache directory. Existing
cache files are rehashed before use; mismatches fail rather than overwrite.

Agarwal et al. (2025), DOI
[10.1038/s41586-024-08430-9](https://doi.org/10.1038/s41586-024-08430-9), supplies
the old panel, under the previously recorded source license. Kircher et al.
(2019), DOI [10.1038/s41467-019-11526-w](https://doi.org/10.1038/s41467-019-11526-w),
MaveDB and the author OSF deposit
[10.17605/OSF.IO/75B2M](https://doi.org/10.17605/OSF.IO/75B2M) supply the earlier
promoter sources. MaveDB records declare CC0; the separate OSF raw files are not
relabeled CC0. Their existing producer-use statement and attribution remain in
[the source evidence](../../forensic-sources/README.md), and the new
[acquisition receipt](author-acquisition.json) binds all fifteen exact files.
The fifteen unchanged old author tables are retained under
`raw-author-tables/` (1,002,526 bytes), with both source-service MD5 and SHA-256
checked in [author-md5-verification.json](author-md5-verification.json). The
derived inventory emits DNA identity facts and source keys without numeric
response/count fields; retained raw tables still contain their original old-study
columns. Their terms are kept distinct from the CC0 MaveDB snapshots. UCSC and NCBI
reference attribution remains in the source manifests.

The builder and controls are original Helix code under the repository license.
No producer implementation is incorporated. Inventory completeness refers to the
declared explored panel, loci and source-resolved edits; it does not establish
biological independence, validate any new outcome column or replace the frozen
external-study exclusion and evaluation gates.
