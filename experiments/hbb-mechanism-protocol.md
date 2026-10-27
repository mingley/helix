# HBB reference-to-mechanism comparison protocol

## Goal and status

Reconstruct the named HBB transcript and beta-globin protein from GRCh38,
recompute supported single-base consequences without reading ClinVar consequence
annotations, and quantify coverage and disagreements against frozen annotations.

This is a versioned internal analysis specification, not an externally registered
or independent confirmatory study. The existing clinical classifications and
their first benchmark results have already been inspected. They must not be
described as a pristine test set. The new computation is an annotation-consistency
experiment, not a new pathogenicity predictor or a claim to reverse engineer the
whole genome.

## Frozen inputs and units

- `datasets/clinvar-hbb`: all 54 selected SNVs, unmodified source response,
  existing partition assignments and manifest hashes.
- `datasets/hbb-reference`: versioned transcript, protein, chromosome interval,
  exon/CDS coordinates, raw authoritative annotation and source commitments.
- Use one named transcript and assembly. Record strand, coordinate conventions,
  initiator methionine numbering and protein accession.
- The computational unit is a unique genomic SNV. Alternative alleles at one
  locus share context and are not independent biological replicates.

No new source data or changed annotations may silently replace these snapshots.
A new reference, transcript or mapping rule creates a new experiment version,
with reasons and changed commitments recorded.

## Checks before biological interpretation

1. Verify every input commitment and reference allele.
2. Reconstruct the spliced transcript on the correct strand and reconcile it
   against its versioned sequence.
3. Translate the reconstructed coding sequence, handling terminal stop codon
   conventions explicitly, and reconcile the protein against its accession.
4. Test known HbS (`NM_000518.5:c.20A>T`, precursor `p.Glu7Val`) and HbC
   (`c.19G>A`, precursor `p.Glu7Lys`) substitutions as external positive controls.
   Mature beta-chain naming omits the initiator methionine; do not mix mature
   and precursor positions.
5. Test reverse-complement boundaries, exon junctions, splice-proximal sites,
   untranslated regions, reference mismatches and unsupported variants.
6. Enumerate three single-base alternatives over a disclosed finite region.
   Exhaustive SNV enumeration is a computational sensitivity map, not empirical
   saturation-mutagenesis data.

Failure blocks a mechanistic-consistency claim. Correct translation or a matching
codon consequence does not establish effects on abundance, binding,
polymerization, physiology or clinical outcome.

## Primary comparison

The primary descriptive endpoint is the fraction of compatible coding consequence
annotations reproduced for the named transcript. Record the exact mapping from
engine categories to existing `missense`, `synonymous`, `splice`, `stop`, `intron`,
`utr` and `other` groups before producing aggregate numbers. Do not alter the
mapping to improve observed agreement.

Report:

- Total SNVs, compared SNVs, unsupported/context-incompatible SNVs and reasons.
- Coverage: compared / all targets.
- Agreement: matching categories / compared targets.
- Every disagreement with coordinates, transcript, source annotation, computed
  codon/region and a proposed explanation.
- Coding-only and all-region results separately. Splice proximity is a location
  flag, not proof of disrupted splicing.

Source annotations can use different transcripts or consequence priorities.
Disagreement is a question to investigate, not automatically a biological
contradiction or software defect.

The working null is that the independently recomputed category and source
category are consistent under identical assembly, transcript and grouping
semantics. No random null agreement rate is assumed: arbitrary chance agreement
is inappropriate for selected genes with unequal category frequencies.

## Secondary descriptive statistics

`genome-research-stats` reports nominal 95% Wilson score intervals with the
independent Bernoulli assumption attached. The selected fixture and correlated
alleles do not establish this assumption; intervals describe sample-size
sensitivity, not guaranteed population coverage or the probability of biological
truth.

If revisiting the previous consequence model versus prevalence baseline:

- Reuse the exact frozen predictions and target identities.
- Disclose the analysis as post hoc.
- Report all four paired correctness cells and the exact conditional two-sided
  McNemar/binomial calculation on discordant pairs.
- Null: left-only and right-only correctness are equally likely among
  discordant pairs, conditional on their number.
- Report p-value and accuracy difference without turning a threshold into
  discovery. Selected variants, shared studies and prior inspection prevent
  confirmatory population inference.
- Preserve the logarithmic p-value if probability underflows numerically.

The existing held-out set has 21 benign and 4 pathogenic assertions. Predicting
everything benign gives 84% raw accuracy and 50% balanced accuracy. This
test-set majority fraction is a descriptive caution, not a baseline fitted from
held-out labels. Keep the actual baseline fitted on training records. Report
sensitivity, specificity, balanced accuracy, calibration and individual errors.

## Acceptance and follow-up

Software acceptance requires reference reconstruction, explicit accounting of all
targets, repeated stable outputs and real negative controls. There is no target
agreement percentage to optimize toward. Unexplained discrepancies remain open.

A mechanistic hypothesis needs a discriminating independent measurement:
transcript/splicing data for a splice mechanism, protein abundance for a dosage
mechanism, or a controlled functional assay for an amino-acid mechanism. Clinical
curation and annotation recomputation are not two independent experiments.

A future prediction study must freeze fresh independent experimental units,
endpoint, exclusions, assay units, baselines and model-selection protocol before
accessing outcomes. Agreement in this selected HBB fixture does not establish a
general explanation of the human genome.
