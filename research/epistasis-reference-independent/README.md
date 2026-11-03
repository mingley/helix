# Independent S18 reference-mapping check

This outcome-free check reconstructs the companion
[reference audit](../epistasis-reference-audit/README.md) with separate Python
standard-library readers. It parses the original S18 XLSX as ZIP/XML, identifies
the four-diplotype sequence squares by nucleotide-column patterns, parses the
hg19-to-hg38 chain into aligned blocks, and decodes required intervals of the
hg38 2bit reference directly. It does not import the primary mapping code or its
2bit reader. It reads no S20 workbook or response values.

For all **35,616 source rows**, it reproduces all **7,542 accepted design
groups** and all **1,362 rejected group identities and row sets** from the
design report. For every accepted group it independently compares the complete
reference-ledger object: source rows and edit positions, both implied hg19
starts, every distinct whole-interval chain candidate, the source-sequence hash,
all hg38 mismatch positions and bases, and status. The result is **7,480 valid,
58 full-reference mismatches, and 4 with no single-block map**, with every
ledger field matching. The [receipt](receipt.json) binds the source, reference,
chain, design report, ledger and verification implementation by SHA-256.

Reproduce from the repository root after obtaining the hash-pinned hg38 2bit
reference via `scripts/setup_reference.sh`:

```bash
source scripts/env.sh
python3 -m unittest discover -s research/epistasis-reference-independent -p 'test_*.py'
python3 research/epistasis-reference-independent/verify.py \
  --reference "$HELIX_TOOLS_DIR/reference/hg38.2bit" \
  --receipt /tmp/epistasis-reference-independent-receipt.json
cmp research/epistasis-reference-independent/receipt.json /tmp/epistasis-reference-independent-receipt.json
```

The verifier requires explicit receipt output and raises errors on source hash,
design, or ledger discrepancies even under Python optimization. Synthetic tests
cover both chain strands and byte orders, interval boundaries and gaps, duplicate
and ambiguous mappings, N blocks, wrong diplotype labels, every mapping status,
and an altered full-reference base.

The source, chain, assembly reference, and declared eligibility criteria are
shared with the primary audit. The check does not independently classify each
rejected design's detailed reason, prove source assembly documentation, test
alignment outside the selected chain blocks, or establish outcome availability,
replicate support, endogenous function, or biological interaction. Exact
agreement is computational verification of these stated checks.

## Direct hg19 check

The separate [hg19 audit](../epistasis-hg19-audit/README.md) compares the same
S18 reference oligos directly to the genomic reference on the assembly stated by
the source. `verify_hg19.py` reuses this directory's independent XLSX and 2bit
readers and matches **every field of all 7,542 hg19 records**. All 7,542 match
the pinned hg19 sequence exactly. A key-by-key join with the separately checked
hg38 ledger shows that these include all **58 hg38 sequence-mismatch windows**
and all **4 no-map windows**, as well as the 7,480 hg38 exact matches. These
results do not determine why hg38 mapping or sequence checks failed.

The separate [hg19 receipt](hg19-receipt.json) binds both verification scripts,
the source, design report, hg19 reference, hg19 ledger, previously verified hg38
ledger and its receipt. The hg19 verifier checks that this hg38 ledger is bound
to the independent receipt; reproducing the hg38 check first supplies that
evidence. Synthetic controls also exercise mismatches at either edited site and
in the background; the actual hg19 data produced no such mismatch records.

After obtaining the hash-pinned reference described in the hg19 audit, run:

```bash
source scripts/env.sh
python3 -m unittest discover -s research/epistasis-reference-independent -p 'test_*.py'
python3 research/epistasis-reference-independent/verify_hg19.py \
  --reference "$HELIX_TOOLS_DIR/reference/hg19.2bit" \
  --receipt /tmp/epistasis-hg19-independent-receipt.json
cmp research/epistasis-reference-independent/hg19-receipt.json /tmp/epistasis-hg19-independent-receipt.json
```

The direct hg19 comparison and cross-assembly join are source/design checks,
not measurement availability, replicate evidence, a selective rescue of a
failed hg38 gate, or a biological interaction result.
