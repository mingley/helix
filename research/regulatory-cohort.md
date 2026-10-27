# A metadata-selected regulatory cohort beyond HBB

## Goal and scientific unit

The new question is whether reference/alternate sequence context predicts a
published reporter-expression coefficient on promoter elements withheld from
training, and how much identical-allele coefficients change across measured
libraries or contexts. The endpoint remains a **log2 additive regression
coefficient from a plasmid MPRA**, not clinical pathogenicity or endogenous gene
expression.

The cohort comes entirely from Kircher et al. (2019),
[DOI 10.1038/s41467-019-11526-w](https://doi.org/10.1038/s41467-019-11526-w),
[PMID 31395865](https://pubmed.ncbi.nlm.nih.gov/31395865/). Different promoter loci
can serve as cross-element targets within this study. They are not independent
studies, a random sample of the human genome, or independent validation of
clinical classifications.

## Selection was frozen before new outcomes

A public MaveDB metadata-only search by PMID returned 29 score sets. Selection
used species, promoter title, DNA target, declared endpoint, permission, genomic
coordinates and exact reference matching. No coefficient distribution or model
performance was used for inclusion. The complete discovery response, including
excluded enhancers and other contexts, is retained.

All nine eligible distinct non-HBB human promoter elements were retained, below
the cap of ten. The lexicographically first full accession for a repeated locus
is the primary representative. Five other records were frozen as validation-only
library/time/cell/siRNA contexts; HBB is a previously exposed reference only.

The exact roles and analysis rules were frozen in
[regulatory-prediction-v1.json](../experiments/regulatory-prediction-v1.json),
SHA-256 `c8e97fabef2afe8ccaa5a01a7ef7bb952f9c9f15857e5c540af112a66e78cb58`.
A gate containing that protocol hash, the frozen metadata-manifest hash and all
accession roles was written **before the first new score request**. The gate and
the exact protocol bytes remain in the cohort. This establishes a recorded
internal acquisition sequence, not external preregistration or authentication.
The old HBB scores had already been inspected and are identified accordingly.

## Primary elements

Intervals below are zero-based, half-open GRCh38 coordinates. Every deposited
target matches the plus-reference strand exactly within its bounded downloaded
window; gene transcription strand is a separate quantity.

| Element | Primary MaveDB accession suffix | Cell line | Target bp | Interval | Raw rows | Strict SNVs |
| --- | --- | --- | ---: | --- | ---: | ---: |
| F9 | `00000015-a-1` | HepG2 | 303 | chrX:[139530462,139530765) | 984 | 904 |
| FOXE1 | `00000016-a-1` | HeLa | 600 | chr9:[97853254,97853854) | 2,048 | 1,800 |
| GP1BB | `00000017-a-1` | HEL 92.1.7 | 385 | chr22:[19723265,19723650) | 1,268 | 1,154 |
| HBG1 | `00000019-a-1` | HEL 92.1.7 | 274 | chr11:[5249804,5250078) | 907 | 822 |
| HNF4A | `00000020-a-1` | HEK293T | 285 | chr20:[44355519,44355804) | 977 | 855 |
| LDLR | `00000023-a-1` | HepG2 | 318 | chr19:[11089230,11089548) | 1,083 | 954 |
| MSMB | `00000024-a-1` | HEK293T | 591 | chr10:[46046243,46046834) | 2,019 | 1,773 |
| PKLR | `00000027-a-1` | K562, 24h | 470 | chr1:[155301394,155301864) | 1,776 | 1,407 |
| TERT | `00000031-a-1` | HEK293T | 259 | chr5:[1294988,1295247) | 975 | 777 |

The nine primary loci do not overlap. Duplicate or reverse-complement-identical
whole targets are rejected. This does not prove statistical independence between
promoters or between their shorter shared sequence features.

## Validation-only and previously exposed data

| Accession suffix | Role | Context relationship |
| --- | --- | --- |
| `00000023-a-2` | validation-only | LDLR second source library/biological repeat |
| `00000027-b-1` | validation-only | PKLR 48h versus 24h |
| `00000031-b-1` | validation-only | TERT SF7996 glioblastoma versus HEK293T |
| `00000031-c-1` | validation-only | TERT SF7996 GABPA-siRNA knockdown |
| `00000031-d-1` | validation-only | TERT SF7996 scrambled-siRNA control |
| `00000018-a-1` | previously exposed reference | Existing HBB/HEL 92.1.7 reporter data |

Exact genomic allele identity links repeated measurements. Records at the same
locus must not enter training and be counted as independent held-out elements.
The paired contexts can quantify observed coefficient differences and constrain
what a predictor that ignores context can fit. They do not, by themselves,
identify the biological cause of every difference or establish independent
laboratory replication.

Cell names are normalized only from the source descriptions. Library identity,
time, siRNA and treatment remain separate assay contexts; for example, two LDLR
libraries are not mislabeled as different cell lines. All raw descriptions remain
available.

## Complete source accounting and missingness

The retained 15 score sets contain **18,468 raw rows**:

- 15,699 unambiguous reference-validated SNVs.
- 2,769 literal source equality descriptions.
- All source-provided scores in this particular snapshot are finite.

The nine primary elements contribute 10,446 SNVs; validation-only sets contribute
4,692; the previously exposed HBB reference contributes 561. The source tables
remain byte-for-byte unchanged.

Only strict positive one-based target notation `n.POSREF>ALT`, with unequal A/C/G/T
alleles and matching reference base, becomes an SNV record. Target coordinates
are not reinterpreted as cDNA or chromosome coordinates. Source equalities,
complex alleles and other representations remain explicit source rows with
exclusion reasons. Equality notation is not assumed neutral and is not silently
changed into a deletion.

Representation handling is separate from numeric quality. Non-numeric, missing,
nonfinite or out-of-domain scores/counts have explicit states. Ingestion retains
reference-valid SNVs even if a score or count is unavailable. The frozen analysis
protocol, not acquisition, applies minimum tags, positive read counts and any
other declared quality rule. A row is never silently imputed or dropped.

F9 has 904 of 909 theoretically possible SNVs, GP1BB 1,154 of 1,155, and each PKLR
timepoint 1,407 of 1,410. All target positions appear, but some alternate alleles
are absent. The final pre-acquisition protocol explicitly retains quality-passing
SNVs at incomplete positions rather than deleting an entire position or inventing
its missing alternate.

## Provenance, permission and reproducibility

Each retained metadata response and score CSV has its exact source URL, role,
retrieval time, byte count and SHA-256. MaveDB creation/modification dates remain
inside the raw metadata. A small UCSC hg38 JSON sequence response is retained for
each distinct target window; no whole-genome sequence download was required.

Every score set declares a CC0 1.0 **data** grant. The previously reviewed MaveDB
per-score-set license-scope evidence is retained, and the original producer is
attributed here. Article or platform-software licenses are not substituted for
data permission. Data keep their source terms rather than being relicensed under
the repository's MIT code license.

Normal loading is offline. `genome-regulatory-data::RegulatoryCohort::open` checks
input commitments, metadata-first gate, frozen roles, data license, species,
endpoint, cell context, exact target/reference identity and CSV source accounting.
The Rust and Python loaders independently check row normalization.

```bash
source scripts/env.sh
cargo test -p genome-regulatory-data --locked
python3 scripts/fetch_regulatory_cohort.py verify datasets/regulatory-cohort
```

Refreshing is deliberately a two-stage workflow into a new directory:

```bash
python3 scripts/fetch_regulatory_cohort.py metadata datasets/regulatory-cohort-new
# Review/freeze metadata-derived protocol before requesting new outcomes.
python3 scripts/fetch_regulatory_cohort.py acquire datasets/regulatory-cohort-new \
  --protocol experiments/regulatory-prediction-v1.json \
  --expected-protocol-sha256 c8e97fabef2afe8ccaa5a01a7ef7bb952f9c9f15857e5c540af112a66e78cb58
```

Do not reuse an existing complete snapshot or change its role assignments after
viewing outcomes. Preserve incomplete attempts if acquisition fails; a complete
snapshot is identified by its final manifest. New source versions or changed
inclusion rules require a new experiment version.

## Inference boundaries

These are source-fitted regression coefficients from random-mutagenesis libraries,
which can contain multi-change constructs. They are not isolated endogenous
single-variant measurements. The source's p-value limitations, additive-model
assumptions, barcode/library dependence and absence of per-culture coefficients
remain relevant. No new significance threshold is assigned during ingestion.

Cross-element evaluation also changes cell context for several targets. Shared
study processing, disease-associated element selection and unmodeled reporter
context limit transport to new studies and the rest of the genome. A successful
sequence association would motivate an independently designed perturbation test;
it would not establish that genome function has been reverse engineered.
