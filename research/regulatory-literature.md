# Regulatory sequence priors and novelty boundary

## Purpose and study status

The goal is to test whether a small, interpretable set of regulatory-sequence
features predicts the **reporter-expression effects** in the frozen assay
cohort. It is not to label every motif match as occupied, infer a clinical effect,
or treat an expression association as a new biological mechanism.

The HBB assay was already explored before this round. The four motif identities,
their versions and the scoring rule were selected from prior biology before
this subagent inspected positional MPRA score patterns. They were then frozen
for [the study protocol](../experiments/regulatory-prediction-v1.json).
This is a transparent hypothesis-driven reanalysis with additional evaluation,
not a claim that the entire research history was blind.

## 1. Fixed, outcome-independent motif selection

The source is the official [JASPAR API](https://jaspar.elixir.no/api/v1/).
[The retained motif bundle](../datasets/regulatory-motifs/README.md) contains raw
responses, exact normalized counts, source hashes, license evidence and the
scoring contract.

| Profile | Biological reason for selecting it | Provenance that must not be overstated |
| --- | --- | --- |
| GATA1 MA0035.4 | An erythroid transcription-factor prior, included to allow a motif hypothesis rather than assume its binding. | JASPAR records human taxid 9606 and ChIP-seq. This does not establish expression/occupancy in the assay cell line or direct regulation of every target. |
| KLF1 MA0493.3 | EKLF/KLF1 and beta-globin CACCC regulation are established prior biology. | The frozen matrix is annotated human/ChIP-seq. The original EKLF discovery cited below was murine; these are distinct provenance claims. |
| NFYA MA0060.3 | CCAAT boxes and NF-Y are established globin/promoter priors. | NF-Y is a heteromeric complex. “NFYA motif” is the database label, not evidence that isolated NFYA produces the assay effect. |
| TBP MA0108.3 | TATA/core-promoter sequence is a prior independent of observed assay outcomes. | JASPAR lists vertebrates but no species or experimental type. The source publication analyzes eukaryotic promoter sequences. It is not labeled a human-only ChIP/affinity experiment. |

The versioned accessions are the actual selection. A citation to JASPAR 2024
does not imply these profiles all belong exclusively to one named release or
are the latest available versions. The official About page explicitly licenses
JASPAR under CC BY4.0; the source page and API files are preserved unchanged.

### Prior evidence already establishes much of the expected pattern

- Miller and Bieker identified erythroid EKLF and showed binding/activation at
  a beta-globin CACCC element in their studied murine systems. This is established
  biology, not a discovery produced by a new feature vector. [R02]
- Marini et al. investigated both proximal and distal beta-globin CACCC boxes,
  reporting differing binding/activation behavior and a contribution of both
  elements to maximal EKLF-dependent stimulation in the studied constructs.
  A new association near these elements must be compared with that history.
  [R03]
- NF-Y/CCAAT regulation of globin promoters has a substantial literature.
  Developmental and promoter context matter; a motif effect need not transfer
  with one universal coefficient between HBB, HBG and unrelated genes. [R04, R05]
- The TBP profile's original reference derived matrices from 502 eukaryotic
  promoter sequences. Its matrix is a useful sequence prior, not a direct
  biochemical measurement of each HBB allele. [R06]

These priors do not guarantee that four profiles are sufficient for nine other
promoters, other cell types or alternate assay architectures.

## 2. Exactly what the motif feature means

The frozen transformation is:

```text
P(base,column) = (count(base,column) + 0.25) /
                 (sum_base count(base,column) + 1.0)
weight(base,column) = log2(P(base,column)/0.25)
window_score = sum_column weight(base_at_column,column)
delta = max_alternate_window_score - max_reference_window_score
```

The total pseudocount is one per column, apportioned over the uniform background;
it is not one per base. All complete motif-width windows overlapping the variant
are scanned in both orientations. Both alleles use the same candidate
coordinates/orientations, while their maximizing windows can differ.
The protocol fits four signed deltas only; reference maxima may be descriptive,
but are not silently added as fitted covariates. No motif-score threshold or
outcome-selected radius is used.

The widest motif is 11 bases, so an SNV feature can use the union of reference
bases within ten positions on either side. The guard of 20 requires centers
more than 20 bases apart for strict train/test sequence-support separation.
This sequence guard is not study/cohort independence.

A delta is a change in **motif log-odds similarity**, measured in bits. It is not
a change in binding free energy, TF occupancy, transcription, or clinical risk.
Independent-column PFMs omit factor abundance, cooperative binding, spacing,
DNA shape, chromatin, methylation and many sequence interactions. ChIP-derived
profiles can also reflect source-context occupancy rather than an isolated
TF's intrinsic affinity.

## 3. What Kircher 2019 already reported

The study from which the assay cohort comes already performed saturation
mutagenesis of regulatory elements and compared measured effects with
computational annotations, including JASPAR motif predictions. It reported:
[R01]

- Incomplete overlap between experimentally supported binding annotations and
  motif predictions.
- Difficulty modeling gain-of-binding effects from scans of the reference
  sequence alone; TERT activating mutations were a concrete example.
- Large differences in predictive performance across regulatory elements,
  including poor average agreement and context-dependent gains.
- An orientation-dependent region near the minimal promoter in the SORT1
  reporter architecture, with a candidate EBF1 motif interpretation.
- Better directional information from allele-aware sequence models in some
  elements than from reference-position annotations.
- An explicit TERT/GABPA analysis of reference/alternate PWM score changes,
  and TERT cell-dependent effects around a candidate E2F repressor region.

Therefore, none of the following would by itself be a new biological discovery:
“motifs matter,” “creating a motif can affect expression,” “prediction differs
between promoters,” “reference-only motif scans miss gains,” or “reporter context
can change an effect.” Those are prior published observations.

Allele-aware PWM disruption is also established computational methodology:
*motifbreakR* (2015) compares allele-specific motif information and *atSNP*
(2015) assesses variant changes in motif-match scores under explicit statistical
background models. Our deliberately simpler uniform-background delta is not a
new general method for discovering TF binding. [R08, R09]

TERT gain-of-binding is especially unsuitable as a novelty claim: direct GABP
recruitment to recurrent mutant TERT promoters was reported in 2015, before the
assay paper. The fixed four-profile set does not even include ETS/GABP motifs;
a failure in that context is not evidence that the published mechanism is
mysterious or false. [R10]

Our narrower contribution can be an auditable, independently reproduced
quantification under a particular frozen scoring rule, feature set, split and
context test. It may produce a useful new *computational observation*, but that
must be stated at the level actually supported by its controls and uncertainty.

## 4. Reporter architecture: what is resolved and what remains unknown

### Different coordinate and physical orientations

The MaveDB HBB target is represented on the GRCh38 forward strand. Native HBB
transcription is on the minus strand. The genomic coordinate representation
alone does not tell us how the promoter fragment was cloned.

The reported design, supplement and deposited backbone sequences provide a
stronger evidence chain:

1. Supplementary Table 1 identifies the HBB MPRA backbone as pGL4.11c and the
   cell line as HEL92.1.7.
2. Supplementary Table 18 gives directional HBB amplification oligos. The
   forward oligo's target-specific suffix matches the start of the **reverse
   complement** of the MaveDB target, while the reverse oligo's suffix matches
   the start of the plus-strand target.
3. The MCS ordering in deposited pGL4.11b/MK484103.1 and
   pGL4.11c/MK484104.1 is consistent with the native-minus-orientation HBB
   fragment upstream of luciferase.
4. The main text reports that the reporter is inserted between the regulatory
   element and the random barcode; the barcode then lies in the reporter
   3′ UTR. [R01]

This resolves the **reported DNA design orientation**; it is not a new
sequencing verification of every individual library molecule. The evidence and
source hashes are in
[construct-review.json](../datasets/regulatory-motifs/construct-review.json).

Both-strand DNA motif scanning is invariant to reverse-complementing the entire
reference/alternate target pair. It avoids mistaking the reference-coordinate
strand for the biologically meaningful orientation of one occupied site.

### Native 5′ UTR inclusion is not a mapped reporter transcript

The insert covers the native 50-base HBB 5′ UTR as well as upstream promoter
sequence. Its physical inclusion in the DNA design does **not** establish the
reporter transcript's exact start or complete RNA sequence. Native-TSS initiation
would predict transcription of that portion, but the assay's barcode RNA/DNA
readout does not map 5′ ends, cap sites, splice isoforms or allele-specific
initiation choices.

No uAUG, uORF, RNA-stability or translation feature is justified merely by
intersecting the target with native RefSeq annotation. The reporter endpoint is
RNA abundance relative to DNA-derived barcode counts, not a direct measurement
of luciferase protein, beta-globin protein or hemoglobin physiology. An RNA-based
mechanism would require the exact final construct and appropriate transcript
mapping/measurements. The frozen feature protocol therefore adds no such RNA
feature.

## 5. How eventual findings will be reviewed for novelty

A candidate needs a precise claim before a novelty search: sequence/allele,
target and cell context, scored endpoint, direction of contrast, comparison
model, uncertainty and the controls that survived. We will then distinguish:

| Category | Appropriate interpretation |
| --- | --- |
| Established mechanism recovered | Reproduction of known promoter/TF biology, with a new implementation or dataset-specific estimate |
| New computational reanalysis | A new quantified property of this frozen analysis, such as transfer failure after explicit controls; not automatically a new mechanism |
| Unvalidated novelty candidate | A specific pattern not located in the bounded literature search and not explained by current controls; requires independent confirmation |
| Negative or inconclusive result | A feature did not improve the selected endpoint under the frozen protocol, or evidence is too uncertain to decide |

Searches will use factor aliases (KLF1/EKLF, NF-Y/NFYA), gene aliases
(HBB/beta-globin), coordinate/allele aliases and the actual reporter/context.
The original assay paper and its supplement are part of the prior art, not an
independent validation source for a claim inferred from those same scores.

“No match found in this search” is not a proof of novelty. A small motif
association does not establish TF binding; a new regression coefficient does not
establish causality; and repeated holdouts from one original experiment do not
become independent biological replication.

The final candidate-specific assessment belongs in
`research/regulatory-novelty-review.md` after the frozen analyses finish.

## Sources and access scope

[R01] Kircher et al. *Saturation mutagenesis of twenty disease-associated
regulatory elements at single base-pair resolution.* Nature Communications
(2019). [PMID 31395865](https://pubmed.ncbi.nlm.nih.gov/31395865/),
[doi:10.1038/s41467-019-11526-w](https://doi.org/10.1038/s41467-019-11526-w).
Main text and supplementary PDF were read through the legitimate Europe PMC
archive. The article is CC BY4.0. Article rights and a separately deposited
score-set license are distinct.

[R02] Miller IJ, Bieker JJ. *A novel, erythroid cell-specific murine transcription
factor that binds to the CACCC element and is related to the Krüppel family of
nuclear proteins.* Molecular and Cellular Biology (1993).
[PMID 7682653](https://pubmed.ncbi.nlm.nih.gov/7682653/),
[doi:10.1128/MCB.13.5.2776-2786.1993](https://doi.org/10.1128/MCB.13.5.2776-2786.1993).
Bibliographic record and abstract read; murine study context retained.

[R03] Marini et al. *The distal beta-globin CACCC box is required for maximal
stimulation of the beta-globin gene by EKLF.* British Journal of Haematology
(2004). [PMID 15384985](https://pubmed.ncbi.nlm.nih.gov/15384985/),
[doi:10.1111/j.1365-2141.2004.05153.x](https://doi.org/10.1111/j.1365-2141.2004.05153.x).
Bibliographic record and abstract read.

[R04] Martyn GE, Quinlan KGR, Crossley M. *The regulation of human globin promoters
by CCAAT box elements and the recruitment of NF-Y.* BBA Gene Regulatory
Mechanisms (2017). [PMID 27718361](https://pubmed.ncbi.nlm.nih.gov/27718361/),
[doi:10.1016/j.bbagrm.2016.10.002](https://doi.org/10.1016/j.bbagrm.2016.10.002).
Review abstract read; used as background, not a new experimental result.

[R05] Fang et al. *Developmentally specific role of the CCAAT box in regulation
of human gamma-globin gene expression.* Journal of Biological Chemistry (2004).
[PMID 14645237](https://pubmed.ncbi.nlm.nih.gov/14645237/),
[doi:10.1074/jbc.M306241200](https://doi.org/10.1074/jbc.M306241200).
Abstract read; developmental/construct limits retained.

[R06] Bucher P. *Weight matrix descriptions of four eukaryotic RNA polymerase II
promoter elements derived from 502 unrelated promoter sequences.* Journal of
Molecular Biology (1990).
[PMID 2329577](https://pubmed.ncbi.nlm.nih.gov/2329577/),
[doi:10.1016/0022-2836(90)90223-9](https://doi.org/10.1016/0022-2836(90)90223-9).
Abstract read; this is the TBP profile's cited source, not proof of a human-only
experimental matrix.

[R07] Rauluseviciute et al. *JASPAR 2024: 20th anniversary of the open-access
database of transcription factor binding profiles.* Nucleic Acids Research
(2024). [PMID 37962376](https://pubmed.ncbi.nlm.nih.gov/37962376/),
[doi:10.1093/nar/gkad1059](https://doi.org/10.1093/nar/gkad1059).
Database abstract and exact profile metadata read; raw profile/license evidence
is retained in the motif bundle.

[R08] Coetzee SG, Coetzee GA, Hazelett DJ. *motifbreakR: an R/Bioconductor
package for predicting variant effects at transcription factor binding sites.*
Bioinformatics (2015).
[PMID 26272984](https://pubmed.ncbi.nlm.nih.gov/26272984/),
[doi:10.1093/bioinformatics/btv470](https://doi.org/10.1093/bioinformatics/btv470).
Abstract read; cited as methodological prior art, not as an implemented dependency.

[R09] Zuo C, Shin S, Keleş S. *atSNP: transcription factor binding affinity
testing for regulatory SNP detection.* Bioinformatics (2015).
[PMID 26092860](https://pubmed.ncbi.nlm.nih.gov/26092860/),
[doi:10.1093/bioinformatics/btv328](https://doi.org/10.1093/bioinformatics/btv328).
Abstract read; its richer statistical/background model is not claimed to be
reproduced by our uniform-background feature.

[R10] Bell et al. *The transcription factor GABP selectively binds and activates
the mutant TERT promoter in cancer.* Science (2015).
[PMID 25977370](https://pubmed.ncbi.nlm.nih.gov/25977370/),
[doi:10.1126/science.aab0015](https://doi.org/10.1126/science.aab0015).
Abstract read; cancer/cell context and cooperative ETS mechanisms are retained.

### Bounded literature-search record

The prior-art pass used Europe PMC title/identifier searches for `JASPAR 2024`,
`CACCC` with `beta-globin`, `CCAAT` with `globin`, the source PMIDs 7682653 and
2329577, `motifbreakR`, `atSNP`, and the exact Bell GABP/TERT title. Kircher's full
text and the HBB design metadata in its supplement were then checked directly.
This is a documented, bounded search rather than a systematic review of all
possible synonyms, unpublished results or the complete regulatory literature.
It does not establish that any later unmatched candidate is globally novel.
