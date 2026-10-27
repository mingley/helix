# What the current HBB evidence can and cannot tell us

This is a **post-hoc descriptive audit** of the already explored demonstration
fixture, not a fresh validation set. It uses the unchanged, hash-pinned NCBI
ESummary response in `datasets/clinvar-hbb/source.json` and the source request,
retrieval time and SHA-256 commitments in `manifest.json`.

The implementation is
`genome_evidence_audit::audit_dataset(&Path) -> Result<AuditReport>`. It first runs
the dataset's complete provenance/reference/partition validation, then rechecks
the raw-source commitment before auditing all listed ESummary records. It does
not refresh ClinVar or silently replace any assertion.

## Selection is doing substantial work

The raw response contains 100 records:

| Aggregate classification | Raw records |
| --- | ---: |
| Uncertain significance | 45 |
| Likely benign | 43 |
| Likely pathogenic | 8 |
| Pathogenic | 3 |
| Benign | 1 |

There are **no aggregate conflicting classifications in this response**. This
does not establish an absence of conflicts across HBB, individual submissions,
or later ClinVar versions. The search selects recent hits using a binary-label
query, then the importer filters actual returned classifications. This is an
ascertainment-selected sample, not a census or unbiased clinical cohort.

The 100-to-54 reduction is exactly:

- 45 uncertain-significance records excluded from the binary task.
- One two-variant haplotype excluded because the importer requires one variant
  per record.
- 54 retained biallelic SNVs.

The excluded haplotype is `VCV003777010.1`,
`NM_000518.5(HBB):c.[70G>T;79G>A]`, aggregate pathogenic. Its components include
p.Val24Phe and p.Glu27Lys. Applying that classification to either component by
flattening the first variation-set entry would be an evidence-assignment error.
The existing importer correctly avoids it.

The binary experiment removes nearly half of the source's unresolved questions.
A high binary score cannot be interpreted as solving those excluded variants.

## Review strength is limited

ClinVar review status describes an assertion-review process. It is not a count
of independent experiments and is not proof of the underlying biology.

| Review status | Raw 100 | Training 29 | Heldout 25 |
| --- | ---: | ---: | ---: |
| Criteria provided, single submitter | 89 | 23 | 23 |
| Criteria provided, multiple submitters, no conflicts | 7 | 5 | 0 |
| No assertion criteria provided | 4 | 1 | 2 |

The conventional ClinVar star mapping assigns these categories one, two, and
zero stars respectively. Neither selected partition contains expert-panel
(three-star) or practice-guideline (four-star) assertions. Unknown future
status strings are reported as `unmapped`, not silently assigned zero.

The raw source lists 113 distinct SCV accessions, none repeated across records.
Distinct accessions identify submissions; one organization or study can supply
many accessions. This ESummary snapshot does not expose the organizational and
experimental-cohort identities needed to establish independence. Multiple
submitters agreeing is not an independent replication result.

Reference: [NCBI ClinVar review status](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/).

## The target is clinically heterogeneous

Training contains 22 likely-benign, one benign and six likely-pathogenic
assertions. Heldout contains 21 likely-benign, two likely-pathogenic and two
pathogenic assertions. Combining likely/definitive categories loses assertion
uncertainty.

Twenty-one training records and 21 held-out records name only unspecified
conditions (`not provided` or `not specified`). Other records name
beta-thalassemia, HBB-related disorder, erythrocytosis, malaria susceptibility
or methemoglobinemia, sometimes together in one aggregate record.

A pathogenicity assertion depends on disease, genotype, inheritance, population
and phenotype. It cannot substitute for a measured effect on splicing, protein
abundance, oxygen affinity, stability or polymerization. Different biological
mechanisms can share one aggregate clinical label.

## Exact loci are disjoint; biological neighborhoods are shared

The locus-hash split prevents identical genomic positions crossing partitions.
Both sets still share HBB, ClinVar curation and a roughly 1.5-kb interval.

- Zero exact cross-partition locus overlaps.
- Minimum cross-partition distance: **3 bp**.
- Two training–heldout record pairs within 5 bp.
- 39 pairs within 50 bp.

The nearest pair is training `VCV003910010.1` at zero-based GRCh38 position
5,226,702 and heldout `VCV003645123.2` at 5,226,705. The next is
`VCV004535916.1` at 5,225,872 and `VCV003910009.1` at 5,225,876.
Nearby variants can share a codon, exon, regulatory neighborhood or study.
Proximity alone does not establish shared participants or experiments.

A subsequent experiment should freeze independent experimental units and,
where appropriate, site/codon/exon groups before selecting a predictor. This
audit does not establish independence of the current submission cohorts.

## Why the demonstration's lone false positive is informative

The already reported model learned from five missense training records: four
pathogenic-class and one benign-class. Its Beta(1,1)-smoothed estimate is
`(4 + 1)/(5 + 2) = 5/7` for every missense target. The same category contains
clinically different variants, so that feature cannot encode every mechanism.

| Field | Frozen disagreement |
| --- | --- |
| Accession | `VCV004685642.1` |
| Transcript description | `NM_000518.5(HBB):c.233A>T (p.His78Leu)` |
| Canonical GRCh38 SPDI | `NC_000011.10:5226658:T:A` |
| Source consequence | Missense variant |
| Clinical assertion | Likely benign |
| Review | Criteria provided, single submitter |
| Condition | Not provided |
| Model probability pathogenic | `5/7` |

Source: [ClinVar variation 4685642](https://www.ncbi.nlm.nih.gov/clinvar/variation/4685642/).
The frozen source dates the assertion's last evaluation to December 3, 2024;
that is assertion metadata, not a newly performed experiment.

This is **not** evidence that the assertion is wrong, that the protein has no
functional effect, or that two measurements contradict one another. It
prioritizes investigation of residue context, abundance, stability, splicing
and disease-specific evidence. Independent codon reconstruction validates an
annotation computation, not a clinical label.

The demo correctly calls 24/25 labels, with 21/25 benign-class targets. The
prevalence baseline calls 21/25 by predicting the common class. This tiny sample
does not establish generalization. Its heldout groups are 18 intron, five
missense, one synonymous and one other. Training has one stop and two UTR
records, but no splice-labelled records. Thus performance does not validate
stop, UTR or splice mechanism prediction comprehensively.

## Questions worth carrying forward

1. Does sequence-derived codon/splice annotation agree with submitted consequence
   annotations? Resolve reference and isoform differences first.
2. Which mechanistic contexts explain different clinical labels in one
   consequence category? Category heterogeneity is not same-variant experimental
   contradiction.
3. How much review strength, disease specificity, study identity and direct
   assay evidence supports each proposed claim?
4. Can a frozen mechanistic predictor improve calibration on a new independent
   assay after endpoint, leakage groups and baseline are fixed?

Models tuned using this audit need new validation data. These already inspected
25 records can no longer be described as untouched evidence for future tuning.
