# HBB promoter reporter-assay snapshot

This is a small, independent **experimental endpoint**, not a validated
pathogenicity predictor or independent replication of every clinical assertion.
It is distinct from the existing ClinVar curation and sequence calculations.

## Source and permission

- MaveDB score set [urn:mavedb:00000018-a-1](https://www.mavedb.org/score-sets/urn:mavedb:00000018-a-1).
- Kircher et al. (2019), “Saturation mutagenesis of twenty disease-associated
  regulatory elements at single base-pair resolution,” Nature Communications 10,
  3583. [DOI:10.1038/s41467-019-11526-w](https://doi.org/10.1038/s41467-019-11526-w),
  [PMID:31395865](https://pubmed.ncbi.nlm.nih.gov/31395865/).
- The actual score-set data are declared
  [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/) by MaveDB.
  `metadata.json` preserves that license object and the producer's data-usage
  policy. Its post-publication requirement to cite the producers is satisfied
  by this attribution.
- `license-scope.json` cites Rubin et al.'s
  [MaveDB platform paper](https://doi.org/10.1186/s13059-019-1845-6), which explains
  that per-score-set Creative Commons licenses apply to downloaded data.
  The paper's own CC BY 4.0 article license is a different grant; it is not
  relabeled CC0.

Raw metadata and the complete score CSV are preserved unchanged. The manifest
records exact retrieval URLs, hashes, byte-level input commitments and the
reference mapping. This data is not relicensed under the code's MIT license.

## Measurement

The study uses error-prone-PCR libraries of regulatory DNA in a plasmid reporter
assay, measured in human **HEL 92.1.7** cells for this score set. RNA and DNA
barcode counts enter a multiple linear regression. The score is the coefficient
for the **log2 variant reporter-expression effect**.

Three cultures transfected with the same constructed library contribute to the
combined model. These are within-study culture/transfection replicates; they are
not three independently constructed HBB libraries or independent studies. Barcode
counts are not replicate counts. The provided CSV does not expose separate
replicate coefficients.

The published analysis requires at least ten tag measurements for a variant.
The primary paper notes that coefficient p-values are statistically inflated
under tag-independence assumptions, correlated counts and unmodeled epistasis.
They are retained as supplied, without a new significance or clinical cutoff.
Random library constructs can contain multiple changes; coefficients estimate
variant effects under the fitted additive model rather than representing every
allele as an individually isolated clone.

## Exact target mapping and representations

The 187-nt target matches the frozen GRCh38 chromosome-11 sequence uniquely on
the **plus** strand at zero-based half-open interval **[5227021, 5227208)**.
The source metadata's start 5227022 / end 5227208 therefore uses one-based
inclusive coordinates. HBB transcription itself is on the **minus** strand.

The fragment contains:

- 50 nt of the named RefSeq transcript's 5-prime UTR: [5227021, 5227071).
- 137 nt upstream of its transcript boundary: [5227071, 5227208).

Consequently a variant within this fragment is not automatically an upstream
promoter variant. A 5-prime-UTR change can affect RNA-level regulation as well as
transcription-associated reporter output. The assay does not by itself separate
these mechanisms.

The raw table has **623** rows:

- **561** unambiguous SNVs: all three alternate bases at each of 187 positions.
- **62** literal equality descriptions, such as `n.1=`.

All source reference bases match the target, accessions and SNVs are unique, and
all 623 supplied scores are finite. Equality rows remain source equality
representations. They are **not** reinterpreted as deletions, collapsed into
wild-type controls, or assumed neutral. Their biological/import-history meaning
needs reconciliation against the original supplementary dataset before use.

## Files and use

| File | Role |
| --- | --- |
| `metadata.json` | Unmodified MaveDB score-set metadata and license |
| `scores.csv` | Unmodified complete CSV response |
| `target.fasta` | Target sequence extracted from metadata |
| `license-scope.json` | Per-score-set permission scope and attribution |
| `manifest.json` | Hashes, sources, mapping, context and validation counts |

Acquisition and missingness rules are specified in
[the experiment protocol](../../experiments/hbb-promoter-assay-protocol.md).
Refresh into a **new** directory:

```bash
python3 scripts/fetch_hbb_promoter_assay.py datasets/hbb-promoter-assay-new \
  --reference datasets/hbb-reference
```

The importer checks permission before downloading scores, exact reference
mapping, source-row count, accession/SNV uniqueness, coordinates, alleles and
score missingness. Verify this retained snapshot without network or writes:

```bash
python3 scripts/fetch_hbb_promoter_assay.py --verify datasets/hbb-promoter-assay \
  --reference datasets/hbb-reference
```

Reporter expression is not endogenous HBB transcription, protein abundance,
oxygen binding, erythrocyte physiology or clinical pathogenicity. The short
plasmid fragment omits distal locus-control and chromatin context. A model trained
on these rows would need position-grouped, context-aware evaluation and
independent perturbation data before stronger conclusions.
