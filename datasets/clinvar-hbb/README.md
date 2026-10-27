# ClinVar HBB demonstration dataset

This fixture contains 54 public, variant-level HBB SNVs on GRCh38 / RefSeq
`NC_000011.10`, not patient genomes or identifiable individual data.

## Provenance

`source.json` is the unmodified NCBI ESummary response for the 100 most recent
search hits returned by the query in `manifest.json`. `reference.fasta` is the
exact NCBI EFetch chromosome-11 interval needed by the selected variants. The
manifest records retrieval time, query, request URLs, reference interval and
SHA-256 hashes for every consumed input.

Versioned VCV accessions, germline classifications and review status remain
available for inspection. The loader checks each selected label and SPDI against
the source snapshot. The importer validates all reference alleles.

Sources:

- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)
- [ClinVar review status](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/)
- [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/)
- [NCBI policies and data use](https://www.ncbi.nlm.nih.gov/home/about/policies/)

Source data are not relicensed under this repository's MIT code license.

## Selection and labels

Retain single-variant records whose canonical SPDI refers to `NC_000011.10` with
unequal unambiguous single-base alleles. Remove duplicate representations.

| Aggregate ClinVar classification | Binary target |
| --- | --- |
| Benign / Likely benign | `benign` |
| Pathogenic / Likely pathogenic | `pathogenic` |

Exclude uncertain, conflicting and all other descriptions. Search hits can contain
these records; filtering uses the actual returned classifications. Combining
"likely" and definitive assertions is a demonstration simplification, not an
evidentiary upgrade. Aggregate labels lose disease, inheritance and zygosity
distinctions.

The sole feature is the existing public molecular-consequence group. The importer
prioritizes splice, stop, missense, synonymous, intron, UTR, then other. These are
known annotations, not a discovered mechanism. Targets are absent from the VCF.

## Frozen locus split

```text
sha256("NC_000011.10:<zero-based position>")[0] % 2
0 -> training
1 -> heldout
```

All alternate alleles at a position stay together. Assignment does not inspect
labels. This snapshot has 29 training and 25 held-out records. The loader rejects
overlap, duplicate variants, missing labels, incorrect partition assignments,
reference mismatches and source inconsistencies.

These are disjoint loci, **not independent experimental cohorts**. They can share
submitting organizations, studies and biases. The agent's input API excludes
held-out annotations; this is not filesystem isolation for tool-using agents.

## Model and interpretation

For `k` pathogenic among `n` training examples in a consequence group:

```text
P(pathogenic) = (k + 1) / (n + 2)
```

This is a Beta(1,1)-smoothed binomial estimate. Unseen groups use training
prevalence and are marked extrapolated. The baseline always uses prevalence.
Both predict pathogenic at probability >= 0.5 on identical targets.

These are not calibrated clinical risk scores. This fixture checks the software
and audit trail; it does not establish biological validity, independent
replication, superiority to standard genomic predictors, or discovery.
