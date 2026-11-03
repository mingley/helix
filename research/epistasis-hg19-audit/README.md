# S18 design windows compared directly with hg19

This **design-only** audit compares every sequence-eligible reference/reference
(`RR`) insert from the frozen [S18 design audit](../epistasis-source-audit/README.md)
against UCSC hg19. It reads no S20 workbook, response measurements, or other
outcome values. Its input cohort is the previously checked 7,542 windows in
2,156 pair identifiers, including the windows that do not pass the separate
whole-insert hg19-to-hg38 comparison.

All **7,542 RR inserts exactly match hg19 across all 200 bases**. This includes
the 58 whose mapped full inserts differ from the hg38 reference and the four
for which the separate comparison finds no whole-insert single-block chain
mapping. For the 58, the observed reference-base differences between the two
assemblies are consistent with their original hg19 design sequence; the result
does not establish why the assemblies differ. See the separate
[hg38 mapping audit](../epistasis-reference-audit/README.md) for its full rules
and per-window records.

## Method and artifacts

The source workbook, the earlier S18 eligibility ledger, the hg19 reference,
and the twoBit reader are bound by SHA-256 before processing. For each eligible
window, the code verifies that the two one-based source coordinates and the
ledger's zero-based insert offsets independently imply the same 200-base,
forward-strand interval. It binds the `RR` row back to its source pair, center,
window and distance, checks the source alleles, and compares every base with
hg19. It does not replace or edit any sequence. A mismatch would be retained in
the per-window ledger with its source and hg19 bases and genomic position.

The committed [`results/windows.jsonl`](results/windows.jsonl) contains every
window, its source `RR` row, interval, source sequence hash, exact-match flag
and any differences. [`results/summary.json`](results/summary.json) records the
counts and input hashes; [`results/manifest.json`](results/manifest.json) binds
the outputs, implementation and reader. Reference coordinates in the ledger
are zero-based, half-open; any mismatch positions are one-based.

An [independent verifier](../epistasis-reference-independent/README.md) parses
the S18 XML and the twoBit encoding separately, reproduces all 7,542 ledger
records, and joins them to the hg38 audit by window key. Its
[receipt](../epistasis-reference-independent/hg19-receipt.json) records all
7,480 hg38-valid, 58 hg38-mismatch and four no-map windows as exact hg19 matches.

This check trusts the separately audited selection of four-construct squares;
it independently rereads only the `RR` row for each selected square. Agreement
does not establish genome-wide uniqueness, outcome availability, common
normalization, replicate support, or a biological interaction. A failed hg38
comparison is not automatically an error in a source designed against hg19.
Same-source computational checks are not independent biological replication.

## Reference provenance and replay

The external reference is UCSC's
[`hg19.2bit`](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit),
816,241,703 bytes. UCSC's
[`md5sum.txt`](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/md5sum.txt)
lists MD5 `bcdbfbe9da62f19bee88b74dabef8cd3`; the downloaded file's SHA-256 is
`797bdb4391cbe23bd30268d41d0d86676fd7b36b7f8e76797b1882d585766fa3`.
The approximately 816 MB file is **external to Git**. Its hash, rather than
the download name or a live service response, binds this analysis. The S18
workbook's publisher and license are documented in the S18 design audit.

From the repository root, with the configured research environment and a new
output path:

```sh
source scripts/research-env.sh
mkdir -p "$HELIX_TOOLS_DIR/reference"
# Keep an existing reference; the audit will reject an incorrect hash.
if [ ! -e "$HELIX_TOOLS_DIR/reference/hg19.2bit" ]; then
  curl --fail --location \
    'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit' \
    --output "$HELIX_TOOLS_DIR/reference/hg19.2bit.download" &&
  printf '%s  %s\n' \
    '797bdb4391cbe23bd30268d41d0d86676fd7b36b7f8e76797b1882d585766fa3' \
    "$HELIX_TOOLS_DIR/reference/hg19.2bit.download" | sha256sum --check - &&
  mv "$HELIX_TOOLS_DIR/reference/hg19.2bit.download" "$HELIX_TOOLS_DIR/reference/hg19.2bit"
fi

python research/epistasis-hg19-audit/audit.py \
  --reference "$HELIX_TOOLS_DIR/reference/hg19.2bit" --out /tmp/s18-hg19-replay
cmp research/epistasis-hg19-audit/results/windows.jsonl /tmp/s18-hg19-replay/windows.jsonl
cmp research/epistasis-hg19-audit/results/summary.json /tmp/s18-hg19-replay/summary.json
cmp research/epistasis-hg19-audit/results/manifest.json /tmp/s18-hg19-replay/manifest.json
python -m unittest discover -s research/epistasis-hg19-audit -p 'test_*.py'
```

The audit refuses an existing output directory and mismatched source or
reference hashes. The small synthetic tests exercise coordinate agreement,
source identity, explicit difference retention and unequal sequence lengths;
they do not replace the full-reference replay.
