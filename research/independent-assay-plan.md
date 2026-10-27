# Independent functional evidence: verified acquisition plan

## Decision

Use the **HBB promoter MPRA from Kircher et al. 2019** as the next same-gene
experimental endpoint. It is available now through MaveDB, has a clear
score-set CC0 grant, fits our frozen reference exactly, and exposes regulatory
mechanisms that codon consequence annotations cannot capture.

Acquisition is complete; it has not been used to fit a predictor, choose a score
threshold, establish predictive superiority or claim discovery. The frozen files
and their representation/missingness accounting are in
[datasets/hbb-promoter-assay](../datasets/hbb-promoter-assay/README.md).
The analysis boundary was written in
[the protocol](../experiments/hbb-promoter-assay-protocol.md) before inspecting
outcome patterns. Initial availability checks had exposed the header and a few
example rows; this is disclosed rather than called an untouched holdout.

## Search evidence and false positives

The current ProteinGym substitution metadata was successfully read from its
[official repository](https://raw.githubusercontent.com/OATML-Markslab/ProteinGym/main/reference_files/DMS_substitutions.csv).
It contains no HBB/hemoglobin assay. This is a catalog-specific observation, not
proof that no HBB functional assay exists.

The exact ProteinGym file revision, byte-level hash and extracted availability
metadata are retained in [the source check](assay-sources/availability.json).

MaveDB's public score-set search for `HBB` returned 22 records. Twenty-one refer
to artificial `2HBB_con` sequences or bacterial ribosomal `rplI` constructs,
including a PDB identifier, rather than human beta-globin. Taxonomy, target
sequence and accession validation prevent treating these name matches as HBB
evidence. The remaining record is the human regulatory assay:

| Field | Verified value |
| --- | --- |
| Score set | `urn:mavedb:00000018-a-1` |
| Study | Kircher et al., Nature Communications 2019 |
| DOI / PMID | `10.1038/s41467-019-11526-w` / `31395865` |
| Cell context | HEL 92.1.7 |
| Endpoint | Additive-model coefficient for log2 reporter-expression effect |
| Target | 187-nt HBB promoter/5-prime-UTR fragment |
| Mapping | GRCh38 chr11 [5227021,5227208), target plus / gene minus strand |
| Raw rows | 623: 561 SNVs and 62 literal equality descriptions |
| Permission | Score-set CC0 1.0; cite the producer's published study |

Exact target mapping was independently reproduced against our 2.5-kb reference,
without looking at score outcomes. The source's one-based inclusive metadata
becomes the zero-based half-open interval above. Protein isoform alignment is not
the relevant mapping for a promoter assay; the transcript boundary is relevant.
The fragment contains 50 nt of 5-prime UTR and 137 nt upstream. Mapping everything
in it to “promoter” would erase a mechanistic distinction.

## What this endpoint can and cannot test

It adds experimentally measured reporter regulation in a named cell context.
It can support later questions about whether a sequence-based regulatory model
orders measured variant-expression effects or identifies region-specific effects.

It does not directly measure endogenous beta-globin expression, heme binding,
polymerization, red-cell survival or clinical disease. A plasmid reporter lacks
distal locus-control and endogenous chromatin context. Random mutagenized
constructs can contain several changes; the fitted additive coefficients do not
establish an epistasis-free mechanism for every individual clone.

The study combines three culture/transfection replicates from the same library.
Those repeat measurements and many barcodes are not independent studies. Source
coefficient p-values have known inflation limitations acknowledged by the paper;
do not turn them into a new biological truth label.

The 62 equality-form rows need comparison with the original supplementary source
before their meaning is used. The importer preserves them rather than guessing
that they are deletions or neutral wild-type controls.

## Proposed next evaluation, not executed

First finish descriptive coordinate/evidence reconciliation. A trained model
requires another versioned plan before using outcomes:

1. Freeze features, model, endpoint, missingness rules and baselines.
2. Group all three substitutions at a position together. Prefer contiguous
   target blocks with a guard interval at least as large as the model's input
   sequence window; random row splits leak local sequence and motif context.
3. With only one 187-nt element, report within-element retrospective performance.
   Do not describe it as transfer to other promoters or the genome.
4. Predefine a training-mean zero-information regression baseline and appropriate
   rank/error metrics. Preserve signs/units rather than selecting a binary cutoff
   against this outcome table.
5. If barcode-level/repeated measurements are acquired later, keep them grouped
   under the same variant. A different barcode from the same library is not a
   new independent validation target.
6. Seek an independent promoter assay, endogenous perturbation or new cell
   context for generalization. Choose units before model selection.

Audit direct overlap with current ClinVar records. Curation can incorporate this
paper, and some predictors may already have trained on these public effects.
Do not use assay-derived clinical assertions or predictor labels as features and
then call the assay an independent validation. Record known training-data
exposure and mark unknown exposure explicitly.

## BRCA1 fallback: real availability, additional caveats

The alternative is also available, but it is not needed to justify the same-gene
HBB endpoint:

- MaveDB `urn:mavedb:00000097-0-2`, “BRCA1 SGE Normalized Scores,” returns a
  complete 3893-row CSV (610,024 bytes) with score, two function-score replicates,
  RNA score and two RNA-score replicates.
- Findlay et al. 2018, [DOI:10.1038/s41586-018-0461-z](https://doi.org/10.1038/s41586-018-0461-z),
  PMID `30209399`, tests 13 RING/BRCT-domain exons (2–5 and 15–23 under the named
  transcript numbering) in the endogenous BRCA1 locus, primarily HAP1-LIG4-KO
  cell viability.
- The function score begins with log2 day-11 variant frequency versus plasmid
  library frequency, corrects editing-position bias using day-5 measurements,
  and normalizes exon/replicate synonymous/nonsense controls. It is not simply
  a day-11/day-5 ratio, a cancer-risk probability or a generic protein-stability
  score.
- The exact source transcript is `NM_007294.3`. Current MaveDB metadata maps it to
  GRCh38, while the original paper also describes hg19-based annotations.
  Reconstruct the original versioned transcript and verify every reference
  allele; do not substitute a current isoform or relabel coordinates without
  liftover validation.
- ProteinGym's verified `BRCA1_HUMAN_Findlay_2018` metadata describes 1837
  single-amino-acid mutants, a 1863-aa reference, growth selection and raw
  `function_score`, using a manual −1 binary cutoff. This protein-level processed
  subset is not the same as 3893 genomic SNVs and can erase RNA/splicing
  differences between nucleotide variants encoding one amino-acid substitution.

MaveDB declares the BRCA1 score set CC0, and the ProteinGym Zenodo release
[10.5281/zenodo.15293562](https://doi.org/10.5281/zenodo.15293562) declares MIT.
However, the original Findlay paper's data-availability statement describes
nonprofit use and a separate commercial license. Do not silently resolve that
historical/current permission conflict on the strength of a catalog label.
Verify depositor authority and the intended reuse before redistributing a
BRCA1 snapshot. No BRCA1 outcome snapshot is added here.

Replicate rows must remain grouped by underlying SNV and locus; exon-blocked or
domain-separated evaluation would be stronger than random variants. The assay
is one study and its processed controls span exons, so replicate/exon separation
does not create a genuinely independent experiment. A fresh laboratory or assay
is needed for that claim.

## References

- [HBB metadata](https://api.mavedb.org/api/v1/score-sets/urn:mavedb:00000018-a-1)
  and [raw scores](https://api.mavedb.org/api/v1/score-sets/urn:mavedb:00000018-a-1/scores).
- [Kircher 2019 primary article](https://doi.org/10.1038/s41467-019-11526-w),
  [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC6687891/).
- [MaveDB score-set license scope](https://doi.org/10.1186/s13059-019-1845-6).
- [BRCA1 normalized-score metadata](https://api.mavedb.org/api/v1/score-sets/urn:mavedb:00000097-0-2).
- [Findlay 2018 primary article](https://doi.org/10.1038/s41586-018-0461-z),
  [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC6181777/).
- [ProteinGym reference metadata](https://github.com/OATML-Markslab/ProteinGym/tree/main/reference_files)
  and [versioned release](https://zenodo.org/records/15293562).
