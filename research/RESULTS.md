# HBB research checkpoint

The [goal](GOAL.md) is a source-grounded partial model of sequence-to-function,
with HBB as a complete bounded example and a primary-human reference atlas as
context. Six subagents developed and independently checked reference data,
mechanisms, biological sources, evidence quality, statistics and genome context.
The regulatory assay was acquired as a second experimental endpoint during that
work. No new biological discovery or whole-genome functional explanation is
claimed.

## What is established computationally

| Check | Observed result | Interpretation |
| --- | --- | --- |
| Reverse-strand reconstruction | Three exons → 628 nt NM_000518.5 | Exact agreement with the frozen transcript |
| Translation | 444 nt CDS → 147 aa NP_000509.1 plus stop | Exact RefSeq and reviewed UniProt sequence agreement |
| Known allele controls | HbS Glu7Val; HbC Glu7Lys | Correct distinct genomic alleles, including initiator Met |
| Existing coarse consequences | 54/54, including 17/17 coding records | Annotation consistency, not clinical truth or measured function |
| Local SNV universe | 4,824 substitutions across 1,608 genomic bases | Every non-reference base substitution in the transcript locus |
| Selected coverage | 53/4,824; one additional upstream ClinVar record | This selected fixture covers about 1.10% of that finite universe |

The exhaustive sequence map contains 957 missense, 325 synonymous, 34 stop-gained,
7 stop-lost, 9 start-lost, 24 canonical splice-core, 552 UTR and 2,916 other
intronic substitutions. These are computed sequence categories, not counts of
validated functional defects. An unchanged canonical motif, intronic location or
synonymous codon does not rule out an RNA/regulatory effect.

An independent Python projection agreed with the Rust engine for seven alleles,
including codons across exon boundaries and two distinct DNA substitutions that
both encode Asn109Lys. See [the crosscheck](reference-crosscheck.json). HbS/HbC and
the highlighted distal/proximal heme-interaction positions are absent from the
selected clinical fixture.

## What changed our interpretation of the clinical benchmark

The 100 retrieved ClinVar records yield 54 selected SNVs after excluding 45
uncertain assertions and one two-variant haplotype. Of the selected records,
42 lack a specific condition. The former heldout contains 23 single-submitter
assertions and two without assertion criteria, with no expert-panel or
multiple-submitter entries.

Loci are disjoint, but training/test sites can be only three bases apart. SCV
identities do not establish independent organizations, studies or cohorts.
The former heldout has now been explored and is development evidence.

The old model's 24/25 versus the prevalence baseline's 21/25 comprises four
model-only correct predictions and one baseline-only correct prediction. The
exact paired p-value is 0.375 under assumptions not established by this selected
fixture; superiority is not established.

Using all already-inspected labels, a deterministic classifier with only the
coarse consequence category can fit at most **51/54** records: missense and UTR
groups contain opposing labels. This is a finite-sample information limit, not
a biological law or future performance estimate. Correct sequence annotations
alone do not supply disease context, protein behavior, RNA fate or penetrance.

## A separately produced regulatory endpoint

The [Kircher2019 / MaveDB HBB assay](../datasets/hbb-promoter-assay/README.md)
provides 561 reference-validated SNV effect estimates across a 187-nt target.
Every position has all three alternate bases. Another 62 literal equality-form
rows remain preserved outside SNV summaries; their meaning is not guessed.

The target exactly matches plus-strand GRCh38 chr11 **[5227021,5227208)**. HBB
transcription is minus-strand. It includes **50 nt of transcribed 5-prime UTR and
137 nt upstream**, so calling every position promoter-only would be incorrect.

Scores are additive-regression coefficients for log2 reporter-expression effects
in HEL 92.1.7 plasmid MPRA. Their descriptive SNV range is −0.57 to +1.19, with
median −0.02. They are not endogenous HBB expression or pathogenicity
probabilities. Random library constructs can contain several changes; the
source documents additive-model and coefficient-p-value limitations. Three
cultures share one library and are not independent studies.

Two exact genomic alleles overlap the clinical fixture:

| Genomic variant, zero-based | Target allele | Source coefficient | Reference context |
| --- | --- | --- | --- |
| chr11:5227063 G>A | n.43G>A | −0.08 log2 | 5-prime UTR; maps to RefSeq c.-43C>T |
| chr11:5227090 C>T | n.70C>T | +0.10 log2 | Upstream of the named transcript |

These joins connect allele identity while preserving distinct endpoints. Clinical
curation may already incorporate the experiment; independent clinical
confirmation is not established. No regulatory predictor was fitted or tuned.
See [the mapping crosscheck](promoter-reference-crosscheck.json) and
[descriptive protocol](../experiments/hbb-promoter-assay-protocol.md).

## Genome context and its limits

The frozen curated RefSeq atlas contains **29,949 primary-chromosome loci**,
representing 29,922 NCBI Gene IDs: 19,427 protein-coding, 9,334 noncoding-RNA and
1,188 pseudogene loci. The 27 extra placements correspond to shared X/Y gene IDs.
Merged transcript spans cover 46.02%, exons 3.01%, and coding exon portions 1.12%
of the primary reference bases. This source-specific annotation omits some gene
classes and alternate contigs. Uncovered sequence is not demonstrated nonfunctional.

HBB's nearby globin-cluster genes and the source-backed BCL11A, HBS1L-MYB, HBA and
HBG relationships show why genotype-to-phenotype reasoning requires regulatory,
developmental and cellular context beyond one coding sequence. Physical proximity
alone does not establish a regulatory relationship.

## Reproduce and challenge the checkpoint

From the repository root, after setup:

```bash
source scripts/env.sh
bash scripts/check.sh
cargo run --locked -- research --out runs/<new-name>
cargo run --locked -- audit runs/<new-name>
node scripts/check_research_viewer.cjs runs/<new-name>/research.html
```

The standard checks pass **55 Rust tests**, strict Clippy, two scoped Lean proof
files and offline raw-source rederivation. The optional Node check covers explorer
interactions, not browser layout. The original demo's audit remains compatible.
Research artifacts bind inputs, protocol, hypotheses, report fields and duplicate
outputs to a revision ledger. Hashes are integrity commitments, not signatures.

The next goal is a versioned regulatory-prediction study with meaningful held-out
experimental units, context-specific endpoints, baselines and uncertainty.
Independent endogenous observations or new experimental/cellular contexts are
needed before broader functional-generalization claims. See
[the five hypotheses](hbb-hypotheses.json),
[mechanistic synthesis](hbb-mechanisms.md) and
[follow-on assay plan](independent-assay-plan.md).
