# What a context-aware sequence model can test

## Recommended question and review status

**Can a small, fixed set of sequence features predict published cell-specific
reporter scores on held-out genomic loci beyond context-specific means?**

**Current limit: the old source columns do not satisfy the published ESS
formula.** All 55,276 admitted records have unweighted three-cell sums above
1e-8, with maximum absolute sum 1.0987267. This exceeds rounding and leaves
their exact preprocessing/biological reference unresolved. The benchmark
therefore concerns the supplied columns labeled ESS; it cannot certify them
as the final paper's stated formula.

The raw-activity proposal, **context-prediction-v1**, is now
[blocked at its semantic source gate](context-v1-source-gate.md). The available
S7 workbook supplies design sequences; S9 supplies three **element specificity
scores (ESS)** rather than the required observed activity means. No v1 fit or
negative predictive result follows from that source failure. A distinct
**context-specificity-v2** targets the actual published ESS columns, with its own
[frozen protocol](../experiments/context-specificity-v2.json) before inspecting
their numeric values or fitting. Source files
were already acquired under v1, and schemas/design metadata were inspected:
this is not an untouched payload or a successful v1.

The [v2 study report](context-specificity-study.md) retains the source failure
alongside a **13.434% primary MSE improvement** on the unchanged columns.
All **243 independent comparison categories pass**, including exact core
features and all 829,140 model predictions. Numerical reproducibility does
not resolve source semantics; the full workflow retains exit 2 for that gate.

This tests context-sensitive sequence dependence in a specified reporter
assay. Success would not identify a causal transcription factor, reproduce
native chromatin, establish a clinical effect, or explain the whole genome.
Failure would limit the selected representation/model under this evaluation,
rather than show that cell context is unimportant.

The initial review used publication text and official metadata before numeric
ESS inspection. Workbooks were acquired under v1; headers/design metadata
preceded v2, while numeric analysis followed its separate freeze. This update
uses post-fit model and source-arithmetic summaries. Published findings,
acquisition, semantic inspection and the observed formula failure are prior
exposure; neither v2 nor the research history is wholly blind. See the
[source inventory](context-literature-sources/README.md) and
[structured evidence notes](context-literature-sources/review.json).

The selected v2 source is the common 60,000-element library in
**Agarwal et al. 2025**, measured in HepG2, K562 and WTC11 human induced
pluripotent stem cells. Same-fragment comparisons avoid changing the DNA when
changing cells. The source paper already models cell specificity; the general
question and multi-task sequence models are prior art. The proposed contribution
is a narrower, auditable evaluation of a fixed small model. Global novelty
is unestablished. [C01]

## The observation includes experimental context

Let `y[i,c]` denote a fragment's **published reporter-activity estimate** in
experimental cell context `c`. Sequence alone is not the complete experimental
input. Cell state, delivery, integration, reporter backbone, promoter, fragment
orientation/boundaries, sampling and processing contribute to what is measured.

For the Agarwal joint library, the documented design includes: [C01, Methods:
“Joint library tested in HepG2, K562, and WTC11 cells”; “Generation of MPRA
libraries”; “Lentiviral infections and DNA and RNA barcode sequencing”]

- A common constructed library of 200-base elements, largely in one orientation:
  sense for promoters and forward for candidate enhancers.
- Lentiviral delivery with three independent infections per experiment, followed
  by collection after three days. Intended multiplicities of infection differ:
  50 in HepG2, versus 10 in K562 and WTC11.
- Candidate elements upstream of a minimal promoter, with RNA-derived barcodes
  measured relative to DNA-derived barcode representation.
- A minimum of ten independent barcodes in source activity processing.
- WTC11 passages 43–49. The authors state cells were not authenticated or
  checked for mycoplasma contamination. These are source limitations, not
  evidence that a particular culture was contaminated or misidentified.

Integration is into the reporter experiment's genomic sites, not restoration
of each fragment's native locus, long-range contacts or endogenous target gene.
Different intended delivery multiplicities reinforce the term **experimental
cell-context contrast**, rather than an isolated causal effect of TF abundance.
RNA/DNA normalization does not prove that every copy-number, integration or
cell-state difference is removed. [C01, C02, C03]

### Distinguish base activity from source specificity scores

Published **base activity** is derived from barcode `log2(RNA/DNA)` values.
Each replicate's activity distribution is normalized to its median, and three
replicates are averaged. These full-library source operations are part of the
endpoint, not training-only transformations implemented here. [C01, Methods:
“Replicates, normalization and RNA/DNA activity scores”]

A difference `y[i,K562] - y[i,HepG2]` therefore compares activity **relative to
each source replicate's library median**. It is not an absolute difference in
transcripts per cell or a direct physical transcription-rate measurement.
A library-wide shift shared by all fragments can be removed by source centering.
Averaging also hides individual replicate estimates.

The paper describes **element specificity score (ESS)** as another
transformation: each cell's joint-library activity values are first
z-score-transformed, then the element's across-cell mean is subtracted. The
final publication identifies
the full ESS table as Supplementary Table 10. Pairwise differences of published
base activities are a different estimand from ESS; neither should silently
replace the other. [C01, Methods: “Calculation of element specificity scores”]

V2 uses the published ESS values unchanged, without reconstructing the source
z-scores or forcing their row sums to zero. The paper's definition is a
per-cell z-score followed by subtraction of the element mean across **all
three cells**, not the mean of the other two. The Methods paragraph does not
specify the SD convention or exact normalization population. A source-wide
construction is documented, but no unverified `ddof`, missing-value treatment
or recovery of raw activity should be invented.

That stated formula implies a zero row sum, two degrees of freedom, and loss
of a common standardized activity component. Under that **conditional theory**,
the shared-DNA mean target is zero, S collapses to G and A to C. It is not the
observed arithmetic of frozen S9: its nonzero sums contradict a rounding-only
explanation. Shared means in unequally scaled contrasts can contain relative
information; they must not be called recovered common activity.

A separate post-hoc arithmetic audit finds rank two, whole-source column
means near zero/sample SDs near one, and a fitted positive weighted constraint:
weights (HepG2, K562, WTC11) approximately (1.024591106, 0.907411107,
1.067997787), with full-source maximum weighted residual 3.055e-14. This is
consistent with additional unequal column scaling after centering, but does
not prove historical preprocessing. Values and fitted targets were not changed.
Both final paper and [preprint](https://pmc.ncbi.nlm.nih.gov/articles/PMC10028905/)
give the unweighted formula; inspected author
loading/prediction code does not resolve its creation history. See the
[study report](context-specificity-study.md) and
[source-semantic review](context-specificity-semantic-review.md) for the
frozen audit and limitations.

Zero cannot be read as inactivity; under the intended contrast endpoint it is
a relative reference, whose exact old-source processing remains unresolved.
Pair-contrast improvement also obeys a balanced-ridge identity and is not an
independent biological corroboration.

Avoiding further outcome normalization, fitting model means/scales on training
data only, and holding out loci can protect our model evaluation. They do not
turn the processed source endpoint into a prospective raw-count pipeline.
This is fixed-source-endpoint prediction. Source processing used the library
including the later evaluation panel; disclose that fact without treating every
published normalization as fitted-model leakage.

### The source panel was selected using earlier outcomes

The joint library is not a random, outcome-naive genome sample. Enhancers were
sampled across activity bins measured in earlier large-scale libraries.
Promoters were selected across bins of mean activity across cells, together
with the 1,000 most variable promoters. The library also contains selected
positive/negative and shuffled controls. [C01, Methods: joint-library design]

This ascertainment can enrich activity ranges and context differences.
Held-out-locus performance describes this selected panel; it does not establish
genome-wide prevalence of cell specificity or performance on arbitrary human
DNA. Preserve source categories and controls for predeclared reporting.
Selecting categories after seeing model performance requires an explicitly
exploratory analysis.

## Falsifiable interpretations

These are mechanistic questions motivating prospective protocols, not results.
QCTX-01's raw-activity contrast requires its stated endpoint, unavailable for
v1; the distinct v2 ESS test must retain its relative-score interpretation.
Predictive
success requires out-of-sample improvement under a declared loss/baseline and
uncertainty appropriate to the generalization unit. A coefficient or attractive
example is insufficient.

| Question | Supporting observation and limited interpretation | Falsifier or alternative | Independent mechanistic validation |
| --- | --- | --- | --- |
| **QCTX-01: Do sequence-by-cell interactions predict same-fragment contrasts?** Compare shared slopes/cell intercepts against a fixed interaction model, holding all contexts of each locus together. | Improved held-out-locus contrast error beyond the shared model and training-mean contrast supports predictable context sensitivity of these features. | No useful improvement; an intercept shift or panel composition explains the gain; overlapping fragments enter train and test. A negative result limits these features/model. | Matched delivery/context perturbations and new measurement of the same sequence set; endogenous tests for transfer outside the reporter. |
| **QCTX-02: Do prespecified TF-family priors add information beyond GC/composition?** Prior candidates include HNF4A/G in HepG2, GATA–TAL1 in K562, POU5F1–SOX2 in WTC11. | Fixed features improve held-out contrast prediction beyond composition, with consistent declared robustness comparisons: a family-level sequence association. | Gain disappears with composition/control accounting, varies across folds, or requires outcome-selected motifs. Several TFs can recognize one motif. | Assayed-cell TF abundance/occupancy, directed motif changes, TF perturbation/rescue, and exclusion of indirect cell-state effects. |
| **QCTX-03: Does a gain depend mainly on source-selected classes?** Preserve promoter, source-cell enhancer, synthetic and shuffled-control metadata. | Improvement remains visible within predeclared biological categories, reporting uncertainty and failed groups. | Controls or activity-enriched categories dominate, restricting the claim to the selected mixture. This need not invalidate measurements. | A new panel sampled independently of prior activity/variance; prespecified sampling weights for population inference. |
| **QCTX-04: Are context contrasts consistent with replicate precision?** Keep separate infection estimates if available and specified. | Predictions and contrast directions are stable across independent infections; uncertainty-aware comparisons retain the gain. | Contrasts are infection-unstable or low-coverage estimates dominate. Barcode counts alone do not estimate all uncertainty. | Additional infections and batch-balanced experiments. Replicate scatter is not automatically a universal biological noise floor. |

For QCTX-01, a schematic model is
`y[i,c] = intercept[c] + shared_features[i] + interaction_features[i,c]`.
The shared term cancels in a within-fragment contrast. This makes the
interaction question more direct than predicting high versus low activity;
it is a statistical decomposition, not a kinetic law. Keep all contexts of a
test fragment outside training. Group or exclude nearby/overlapping fragments,
repeated sequences and reused controls. Three cell labels do not identify an
unrestricted function over all cell types.

QCTX-02's families are established source-study hypotheses, not independently
discovered motifs. A PWM is a sequence prior, not occupancy, binding energy or
transcription rate. Shared recognition, correlated motifs, co-binding and
composition can make individual coefficients non-identifiable. Agarwal et al.
explicitly warn that the expressed TFs examined do not comprehensively cover
family members recognizing each motif. [C01; C03]

## What earlier studies already establish

### Inoue 2017: integrating versus nonintegrating delivery

Inoue et al. designed 2,236 liver-enhancer candidates plus 204 controls with
171-base inserts, selected using HepG2 chromatin/TF annotations. They measured
the same library in HepG2 with integrating wild-type and nonintegrating
mutant-integrase packaging, multiple barcodes and independent infected cultures.
They reported context differences and prediction using sequence/endogenous
biochemical features. [C02]

This shares candidate sequences and a related architecture, but does not
integrate at native enhancer sites. The reported collection schedules and
intended viral copy numbers also differ between conditions. The design does
not isolate “chromatin” as the sole cause. The authors discuss integration-site
effects and use antirepressor elements to mitigate them; mitigation is not
complete-removal proof. [C02, construct sections; Methods: “Lentivirus
packaging, titration, and infection”]

Data accession: **GSE83894**. Metadata lists a 2018 corrigendum,
DOI `10.1101/gr.237321.118`, PMID 29717003. Its text was not successfully
obtained through available official endpoints here. Resolve its scope before
certifying a quantitative reuse; neither assume wholesale invalidity nor ignore
the correction. [C02a]

### Klein 2020: assay, length, orientation and RNA effects

Klein et al. compare nine MPRA strategies using the existing 2,440-element
Inoue array; its 2,236 genomic candidates also underlie orientation and longer
fragment experiments. These are new comparisons with reused sequence selection,
not independent enhancer panels. They vary integration, element/barcode
positions, orientation and extended fragment length. [C03]

The prior art closely matches the proposed question: lasso models use
sequence-derived and biochemical features, and model pairwise assay differences
as residuals after loess adjustment. Features include TF/RNA-binding predictions,
composition and endogenous annotations. “Sequence predicts context-dependent
assay differences” is explicitly published. [C03, feature/regression Methods;
Results: design comparisons]

When a candidate insert is transcribed in a reporter's 3′ UTR, its associations
can reflect RNA stability/splicing as well as transcriptional enhancement.
That rationale does not license RNA/uAUG claims for a different upstream-insert
reporter: actual construct/transcript architecture is required. Correlated
features limit causal interpretation. Orientation independence is an empirical
approximation for the tested enhancer panel; promoter-overlapping fragments
behave differently. [C03]

All nine designs were run simultaneously to reduce batch effects, which does
not erase every delivery/design difference. The paper describes poorer-quality
replicates and exclusions in two mutant-integrase designs. Honor documented
final processing and record deviations. Data accession: **GSE142696**. [C03]

### Alternative multi-cell and factor-perturbation sources

Tewhey et al. test regulatory-variant libraries in lymphoblastoid lines NA12878
and NA19239, and HepG2. Paired reference/alternate sequences sit upstream of a
minimal promoter with barcode readout. This is a useful independent-study
alternative for cell-context **allelic effects**, but eQTL ascertainment and the
variant endpoint differ from joint 200-base fragment activity. Data accession:
**GSE75661**. A 2018 correction is listed as DOI
`10.1016/j.cell.2018.02.021`, PMID 29474912; only correction metadata was
available. Review correction scope and corrected files before reuse. [C04, C04a]

Ulirsch et al. compare K562 with GATA1-overexpressing K562 in a variant MPRA
and follow selected variants with endogenous editing/target-gene studies.
This is a more directed trans-factor contrast than unrelated cells; however,
overexpression can change broader cell state and is not an isolated binding
assay. Erythroid GWAS ascertainment, short allelic fragments, quantile-normalized
processing and the variant endpoint differ from joint-library activity.
Its follow-up illustrates evidence required to move from reporter effects to
endogenous genes and phenotypes. [C05]

### Agarwal 2025 already predicts cell specificity

The common-library source evaluates specificity with multi-task sequence
models, assigns each element across cells to the same cross-validation fold,
and interprets cell-specific motifs. It distinguishes housekeeping/context
families, checks selected TFs against perturbation data, and examines motif
dose/combinations. These concepts are prior art. [C01, Fig. 5; Methods:
“Data pre-processing and model training”; motif interpretation/modeling]

The paper warns that pretrained genome-wide biochemical models may have seen
test-locus annotations during pretraining. A from-scratch small model can make a
clearer exposure claim; this is an evaluation property, not proof of a superior
biological explanation. [C01, “Sequence-based models predict activity”]

## Independence from the original Kircher 2019 benchmark

New measurements, independent sequence selection, unrelated laboratory
provenance and independent software are different independence claims.

| Source | Relationship to Kircher 2019 saturation-MPRA | Boundary |
| --- | --- | --- |
| Kircher 2019 [C06] | Original 20 disease-associated elements; GSE126550 and OSF 75B2M. Earlier HBB/multi-promoter benchmark source. | Another score set here is not a new source-study validation. |
| Inoue 2017 [C02] | Earlier liver-enhancer study; GSE83894; overlapping investigators and assay lineage. | Distinct experiment, but unrelated-laboratory validation and complete sequence/locus disjointness are unestablished. |
| Klein 2020 [C03] | Retests Inoue array; uses pGL4.23c backbone MK484105 deposited by Kircher source. | New design measurements with dependent sequence selection/related construct provenance. Inoue/Klein are not two independent enhancer panels. |
| Agarwal 2025 [C01] | New large-scale/joint measurements, several investigators shared with Kircher. Joint controls explicitly derive from Smith 2013/Ernst 2016; lentiMPRA lineage shared. | Distinct measurement source supported. Unrelated-lab provenance, absence of shared controls and genomic disjointness from the original panel require separate checks. |
| Tewhey / Ulirsch 2016 [C04, C05] | Separate stated author teams/variant ascertainment; cited by Kircher. | Potential independent-study validation subject to correction review, compatible endpoints and exact sequence/locus overlap accounting. Citation alone is not data reuse. |

Before splits, reconcile exact sequence identity, reverse complements,
overlapping intervals, shared loci, controls and repeated measurements. A gene
name alone does not prove identical constructs; a changed record ID does not
prove a new fragment. Preserve full provenance even for excluded controls.

Agarwal explicitly reuses older MPRA controls, including some 171-base controls
extended to 200 bases. Its WTC11 lineage also includes neural-induction controls.
These are stated relationships, not an allegation that the new library
duplicates old measurements. Definitive overlap requires designed sequences
and source records. [C01, library-design Methods]

## Novelty ledger before new outcomes

| Candidate statement | Status |
| --- | --- |
| Integration versus episomal encoding changes activity. | Prior art: Inoue 2017 and Klein 2020. |
| Fragment length/orientation/reporter position can affect measurements. | Prior art, especially Klein 2020; design/panel dependent. |
| Sequence predicts assay differences or cell specificity. | Explicit prior art: Klein 2020 and Agarwal 2025. |
| HNF4-, GATA–TAL1- and POU5F1–SOX2-related features are relevant in these cells. | Established source interpretation. Fixed-profile replication is a representation test. |
| A small fixed model succeeds/fails with strict held-out loci, explicit overlap policy, loss and baselines. | New-to-analysis numerical result on supplied columns, with complete independent numerical reproduction. Source semantics and global novelty remain unresolved. |
| Old S9 violates its stated unweighted ESS formula but obeys a fitted unequal-weight plane. | Reproducible source-arithmetic observation; compatible scaling is not historical proof or a biological discovery. No target repair or refit follows. |
| Gains survive source-category/control/normalization/pretraining-exposure accounting. | Candidate reliability angle requiring frozen choices and data; does not establish new biology. |
| A specific TF directly causes the contrast. | Unsupported by prediction alone; needs discriminating perturbation/occupancy/sequence evidence. |
| This reconstructs endogenous mechanisms, clinical effects or the human genome. | Unsupported by these assays and this review. |

No new biological mechanism is established here. The next goal is authoritative
normalization/raw activity or independently compatible measurements, preserving
the fixed-column benchmark and failed formula check. A weighted-plane fit is
not provenance, and target projection must not force the formula to pass.

## References

- **[C01]** Agarwal V et al. *Massively parallel characterization of
  transcriptional regulatory elements.* Nature (2025).
  DOI [10.1038/s41586-024-08430-9](https://doi.org/10.1038/s41586-024-08430-9);
  PMID [39814889](https://pubmed.ncbi.nlm.nih.gov/39814889/);
  [open article](https://pmc.ncbi.nlm.nih.gov/articles/PMC11903340/).
- **[C02]** Inoue F et al. *A systematic comparison reveals substantial
  differences in chromosomal versus episomal encoding of enhancer activity.*
  Genome Research 27, 38–52 (2017).
  DOI [10.1101/gr.212092.116](https://doi.org/10.1101/gr.212092.116);
  PMID [27831498](https://pubmed.ncbi.nlm.nih.gov/27831498/);
  [open article](https://pmc.ncbi.nlm.nih.gov/articles/PMC5204343/).
- **[C02a]** Inoue F et al. Corrigendum (2018).
  DOI [10.1101/gr.237321.118](https://doi.org/10.1101/gr.237321.118);
  PMID [29717003](https://pubmed.ncbi.nlm.nih.gov/29717003/).
  Correction existence verified; text unavailable here.
- **[C03]** Klein JC et al. *A systematic evaluation of the design and context
  dependencies of massively parallel reporter assays.* Nature Methods 17,
  1083–1091 (2020).
  DOI [10.1038/s41592-020-0965-y](https://doi.org/10.1038/s41592-020-0965-y);
  PMID [33046894](https://pubmed.ncbi.nlm.nih.gov/33046894/);
  [article record](https://pmc.ncbi.nlm.nih.gov/articles/PMC7727316/).
- **[C04]** Tewhey R et al. *Direct Identification of Hundreds of
  Expression-Modulating Variants using a Multiplexed Reporter Assay.*
  Cell 165, 1519–1529 (2016).
  DOI [10.1016/j.cell.2016.04.027](https://doi.org/10.1016/j.cell.2016.04.027);
  PMID [27259153](https://pubmed.ncbi.nlm.nih.gov/27259153/);
  [article record](https://pmc.ncbi.nlm.nih.gov/articles/PMC4957403/).
- **[C04a]** Tewhey R et al. Correction (2018).
  DOI [10.1016/j.cell.2018.02.021](https://doi.org/10.1016/j.cell.2018.02.021);
  PMID [29474912](https://pubmed.ncbi.nlm.nih.gov/29474912/).
  Correction existence verified; text unavailable here.
- **[C05]** Ulirsch JC et al. *Systematic Functional Dissection of Common
  Genetic Variation Affecting Red Blood Cell Traits.* Cell 165, 1530–1545 (2016).
  DOI [10.1016/j.cell.2016.04.048](https://doi.org/10.1016/j.cell.2016.04.048);
  PMID [27259154](https://pubmed.ncbi.nlm.nih.gov/27259154/);
  [article record](https://pmc.ncbi.nlm.nih.gov/articles/PMC4893171/).
- **[C06]** Kircher M et al. *Saturation mutagenesis of twenty disease-associated
  regulatory elements at single base-pair resolution.* Nature Communications
  10, 3583 (2019).
  DOI [10.1038/s41467-019-11526-w](https://doi.org/10.1038/s41467-019-11526-w);
  PMID [31395865](https://pubmed.ncbi.nlm.nih.gov/31395865/);
  [open article](https://pmc.ncbi.nlm.nih.gov/articles/PMC6687891/).
