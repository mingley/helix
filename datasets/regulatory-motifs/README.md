# Fixed regulatory motif priors

This is a small, preselected JASPAR PFM set for interpretable regulatory-sequence
features. It was chosen from promoter/erythroid biology before this subagent
inspected positional MPRA score patterns. The HBB assay has already been explored
at summary level elsewhere in the workspace; this selection does not make the
whole research program pristine or prospectively blinded.

| Accession | Name | Width | Exact source taxonomy / method |
| --- | --- | --- | --- |
| MA0035.4 | GATA1 | 11 | CORE, vertebrates, Homo sapiens taxid9606, ChIP-seq |
| MA0493.3 | KLF1 | 8 | CORE, vertebrates, Homo sapiens taxid9606, ChIP-seq |
| MA0060.3 | NFYA | 11 | CORE, vertebrates, Homo sapiens taxid9606, ChIP-seq |
| MA0108.3 | TBP | 7 | CORE, vertebrates; species and assay type unspecified in the API |

TBP links to human protein UniProt P20226, but a linked protein does not establish
the training species of the matrix. Its primary citation describes comparative
analysis of eukaryotic promoter sequences, rather than a human-only occupancy
experiment. That missing-taxonomy limitation is retained explicitly.

The accessions/versions are fixed. No claim is made that all belong to one named
JASPAR release, that they are currently the latest profiles, or that these factors
are expressed and bound in every assay cell type.

## Data, attribution and rights

Raw API responses are retained unchanged as accession-named JSON files.
`motifs.json` is a metadata/count projection plus an exact scoring contract.
`manifest.json` pins all consumed files, source URLs and retrieval times.
`jaspar-about.html` retains the official JASPAR About page whose footer explicitly
licenses JASPAR under [CC BY4.0](https://creativecommons.org/licenses/by/4.0/).
These source data are attributed to the JASPAR collaborative database and are not
relicensed under the workspace's MIT code license.

Cite the versioned profiles and the database publication:
Rauluseviciute et al., *JASPAR 2024: 20th anniversary of the open-access database
of transcription factor binding profiles*, Nucleic Acids Research (2024),
[doi:10.1093/nar/gkad1059](https://doi.org/10.1093/nar/gkad1059).
This citation documents the database; the explicit frozen accessions identify
the actual profiles used here.

## Outcome-independent scoring

For each base `b` and motif column `j`:

```text
P(b,j) = (count(b,j) + 0.25) / (sum_b count(b,j) + 1.0)
weight(b,j) = log2(P(b,j) / 0.25)
score(window) = sum_j weight(window[j],j)
```

The total pseudocount is one per column, distributed in proportion to a uniform
background. This is not “one pseudocount per base.” Counts already processed by
JASPAR are not reverse engineered into a presumed experimental sample size.

For a reference-matching SNV, enumerate every complete motif-width window inside
the supplied target that overlaps its position. Score both each window and its
reverse complement. Reference and alternate alleles use exactly the same
coordinate/orientation candidate set.

Return the reference maximum, alternate maximum and `alternate_max-reference_max`
in bits. Maxima can occur at different candidate sites; a reported candidate
should retain both locations (and any ties) rather than imply one experimentally
observed binding site. There is no binding threshold, occupancy label, motif
selection by score correlation, learned scan radius or p-value in this feature
definition. A target shorter than the motif produces an explicit missing value.

The widest profile is 11 bp. A feature for a substitution can therefore depend
on reference sequence within ten bases on either side. Strictly disjoint sequence
supports require training/test substitution centers **more than 20 bp apart**.
This does not itself exclude shared assays, cell contexts, libraries or cohorts.

## Interpretation limits

- These are motif log-odds similarity features, not binding energies, measured
  TF occupancy, transcription rates or clinical classifications.
- A ChIP-derived motif reflects its source experiments and can include indirect
  or context-dependent occupancy. NF-Y is a protein complex; the NFYA-labelled
  motif does not imply NFYA alone binds a reporter.
- Independent-column PFMs omit cooperativity, motif spacing, DNA shape,
  nucleosomes, methylation and factor abundance.
- The reference target's genomic strand is not proof of physical reporter insert
  orientation or of the reporter RNA sequence. Both-strand DNA scanning avoids
  asserting a transcription direction; RNA/uAUG mechanisms need separate evidence.
- The availability of the motifs and familiar CACCC/CCAAT/TATA biology is prior
  knowledge. Recovering known promoter sensitivities is a replication or reanalysis,
  not a new biological discovery.

See [regulatory literature](../../research/regulatory-literature.md) for the
biology priors and reporter-architecture assessment.

## Verify or refresh

```bash
python3 scripts/fetch_regulatory_motifs.py --validate-existing datasets/regulatory-motifs
python3 scripts/fetch_regulatory_motifs.py datasets/regulatory-motifs-new
```

The importer checks exact profile IDs, names, widths, collection, explicit human
taxonomy for the first three factors, finite count columns, license evidence,
source agreement and hashes. Refresh writes only a new directory. It never reads
assay outcomes.
