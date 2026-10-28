# Prior work relevant to the four-construct test

This is a focused, non-exhaustive literature note, not a claim of priority.
Regulatory interactions, including human variant interactions, are established
research topics. A result from the next test would address the specific behavior
of this repository's frozen 3-mer-plus-GC model on a specified reporter cohort;
it would not establish the first observation of biological epistasis.

## Closest human variant and reporter precedent

[Siraj et al., *Functional dissection of complex trait variants at
single-nucleotide resolution* (2026)](https://doi.org/10.1038/s41586-026-10121-6)
describe variant-pair experiments using the activities of all four allele
combinations. Their methods estimate interaction effects using DESeq2 and
combine measurements across overlapping 200-bp windows and, in some analyses,
cell types; they explicitly consider covariance between windows. Their methods
also describe recoding alleles according to activity and express interaction as
expected additive effect minus observed double effect. A new analysis must bind
its own reference/alternate orientation and sign and must not infer that a
meta-analysis of shifted windows is the same estimand as a single-window
four-construct contrast. The methods identify Supplementary Tables 20 and 21
for variant-pair analyses. Supplement descriptions and the retained method
extract are in [the source archive](variant-transfer/source/README.md).

Siraj is especially relevant prior art, not an independent new experiment if we
reuse its measurements. Its main single-variant source and some outcomes have
already been inspected in this project: the failed strict v1 read exposed
selected values before stopping, and the separate exploratory v2 study examined
2,153 complete cases. Readers may also have instantiated unselected workbook
values. The published methods and, during this literature review, published
article abstracts and result summaries were read. A future Siraj-based protocol
cannot be described as fully outcome-naive, and exposure to any particular pair
file or result needs its own documented inventory. Prior exposure, the separate
v2 freeze, and missingness limits are recorded in
[the transfer report](variant-transfer-v2-study.md). No new interaction outcome
table was read for this note.

## Motif-combination and other biological precedents

[Agarwal et al., *Massively parallel characterization of transcriptional
regulatory elements* (2025)](https://doi.org/10.1038/s41586-024-08430-9)
analyzed human reporter activity across HepG2, K562 and WTC11 cells, including
motif co-occurrences and modeled super- and sub-multiplicative combinations.
This prior work uses the same underlying research program as our training panel,
so it is neither independent biological validation for that panel nor the same
as a matched factorial contrast of two isolated substitutions in one sequence.

[Adjemout et al. (2025)](https://doi.org/10.1016/j.xgen.2025.100889) report
a human Jurkat-cell luciferase comparison involving an ATP2B4 enhancer haplotype
with five SNPs and an LAX1 promoter SNP, with measurements in triplicate across
three independent experiments. This is an example of studied human regulatory
interaction, but its haplotype-versus-promoter design is not a pair of isolated
SNVs within a single 200-nt insert. The paper also reports a population
association; an assay contrast alone does not establish a population or
clinical consequence.

[Waters et al. (2021)](https://doi.org/10.1186/s13059-021-02574-x) used
quantitative enhancer-FACS sequencing in developing *Drosophila* to compare
wild-type, single-motif and combination-mutant enhancer constructs, and tested
interaction on a log-transformed activity scale. The study illustrates the
four-condition logic and dependence on enhancer context. Its motif mutations can
change several bases or sites in roughly 1-kb fly enhancers; the data therefore
are not a direct test of a human 200-nt single-SNV model.

## Interpretation for the proposed test

For two specified SNVs in one fixed-length construct, the proposed contrast is
`activity(double) - activity(single A) - activity(single B) + activity(reference)`
in the same cell and on a common documented log2 RNA/DNA scale. The repository's
[formal locality result](variant-locality-study.md) establishes an exact zero
for the ideal fixed linear model when the substitutions are at least three bases
apart; it is not an assertion about measured biology. A measured nonzero
contrast in a reporter assay would challenge that model's prediction for that
construct and assay, subject to measurement uncertainty and assay limitations.
A near-zero estimate would not establish absence of interaction without an
appropriately specified precision or equivalence criterion. Neither outcome
alone identifies a transcription factor, endogenous genomic mechanism, or
clinical effect.
