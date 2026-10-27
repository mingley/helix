# Helix

Helix connects typed genomic representations, source provenance, conditional
Lean reasoning and quantitative experiments. The Rust workspace contains
**16 crates** and builds inspectable partial models of human genome function.

The current [regulatory study](research/regulatory-study.md) has a **negative
primary predictive result**: none of five fixed nonconstant models beats the
training-mean baseline across nine held-out promoter elements. It also exposes
source-identity and numerical-reproducibility hazards. No new biological
mechanism or globally first finding is claimed. The next goal is context-aware
prediction on independent data, compared with context-appropriate mean and
sequence-only baselines.

Start with the [research objective](research/GOAL.md),
[study report](research/regulatory-study.md),
[retained artifact snapshot](research/results/regulatory-v1/README.md), and
[novelty review](research/regulatory-novelty-review.md).

## Current evidence

| Result | Scope |
| --- | --- |
| Primary macro MSE: mean 0.11639219; full model 0.11798887 | Nine selected promoter elements from one study; squared log2 reporter-effect units |
| 62 HBB database equality rows reconcile to original deletions | 53 distinct edited targets before QC; 35 rows/35 identities after predefined QC |
| 131,112 Rust/Python predictions agree within 7.73e−12 | Primary losses, fits and controls reproduce; this verifies computation |
| 26 of 240 frozen exact-Spearman checks fail | Full frozen secondary-rank reproduction remains failed |
| Every separately declared prediction grid matches 240/240 checks | Post-hoc numerical sensitivity; changes the rank estimand and does not repair v1 |

The raw deletion identity bound is not a biological noise floor. The motif model's
tiny advantage over chemistry alone is uncertain and still loses to the mean.
Known motif biology, assay-context effects and legacy deletion/equality problems
have substantial prior art.

## Set up and reproduce the study

Validated environment: Linux x86-64, Rust 1.90.0, Lean 4.19.0 and CPython 3.12.
Bash, curl and a C toolchain are prerequisites. No API keys, database, GPU or live
LLM is required. Run from the repository root:

```bash
bash scripts/setup.sh
bash scripts/setup_research.sh
source scripts/research-env.sh
bash scripts/check.sh
bash scripts/check_study.sh runs/study-review
```

Choose a new output prefix. The study command creates separate Rust, Python,
comparison, figures and rank-sensitivity directories. It retains every fixed
model, 99 spatial controls, 2,000 resamples, predictions, fitted preprocessing,
source commitments and failed comparisons.

**The study check returns exit status 2 for the retained secondary-rank
failures**, after producing its figures and diagnostics. Core software/source
checks and successful primary numerical checks are separate from a complete
frozen-study pass. The post-hoc all-grid agreement must not be substituted for
the failed raw exact-rank endpoint.

See the [independent implementation](research/independent-regulatory-replication/README.md),
[frozen protocol](experiments/regulatory-prediction-v1.json),
[statistical review](research/regulatory-statistical-review.md), and
[rank sensitivity](research/rank-sensitivity-review.md).
The published snapshot supports inspection without refitting.

Setup installs tools outside the checkout, under `/workspace/.tools` in cloud
tasks or a local data directory otherwise. Set `HELIX_TOOLS_DIR` before activation
to choose another location. Python research dependencies are hash-pinned in
`requirements/research.lock`. Setup needs network access; normal analysis uses
frozen local inputs. Existing run directories are never overwritten.

## Explore the HBB mechanism map

```bash
source scripts/env.sh
cargo run --locked -- research --out runs/hbb-review
cargo run --locked -- audit runs/hbb-review
```

This reconstructs HBB's named 628-base transcript and 147-residue protein,
enumerates 4,824 SNVs over its transcript locus, independently recomputes
consequences for 54 ClinVar records, audits their evidence, and places the locus
in a primary-human-chromosome annotation context. It integrates the HBB reporter
assay while keeping clinical labels and molecular endpoints separate.

The output contains a readable narrative, JSON annotations and saturation map,
source commitments, a revision ledger, scoped Lean proofs, and a standalone
`research.html` explorer. Open the HTML locally; no server or remote scripts
are needed. An optional Node.js interaction smoke check is:

```bash
node scripts/check_research_viewer.cjs runs/hbb-review/research.html
```

Read the [biological synthesis](research/hbb-mechanisms.md),
[falsifiable questions](research/hbb-hypotheses.json),
[evidence audit](research/clinvar-audit.md),
[genome atlas](research/genome-atlas.md), and
[reporter dataset card](datasets/hbb-promoter-assay/README.md).

The reported reporter DNA design includes native HBB sequence in native-minus
orientation, while the coordinate target is plus-genomic. Its exact reporter
transcript start/isoform remains unmeasured. Native UTR annotation alone cannot
establish an RNA-processing or uAUG mechanism.

## Original clinical demo

```bash
cargo run --locked -- demo --out runs/clinical-review
cargo run --locked -- audit runs/clinical-review
```

This reproduces the original 29/25-record clinical partition. Both partitions
have since been explored and are development data. The earlier 96% versus 84%
accuracy difference does not establish superiority or independent validation;
its exact paired comparison is p=0.375 under assumptions not established here.

Missing Lean is an error. `--skip-proof` is an explicitly unverified demo mode,
and the saved report records that status. Clinical classifications are curated
assertions, not direct measurements of molecular function.

## Scientific contracts

| Layer | Establishes within its scope | Remaining boundary |
| --- | --- | --- |
| Rust contracts | Assembly separation, reference checks and valid representations | Empirical truth |
| Evidence registry | Matching source-backed support and explicit context | Experimental quality and study independence |
| Lean | Conclusions conditional on stated premises and finite computational obligations | Truth of biological premises |
| Statistical evaluation | Performance under the disclosed data, splits and metrics | Population generalization or causal mechanism |
| Independent implementation | Agreement of computations actually checked | Independent laboratory replication |
| Integrity audit | Hash and ledger consistency | Authentication or biological validity |

Proof status never upgrades empirical evidence. Coordinates are zero-based
inter-residue intervals internally and one-based in VCF. Display keys such as
`GRCh38:11:5225592:C>A` are local identifiers; general GA4GH VRS normalization,
computed identifiers, indel support and assembly liftover remain outside the core
variant contracts. The deletion analysis is a source-specific forensic crosswalk.

## Workspace

| Crate | Responsibility |
| --- | --- |
| `genome-types` | Assembly-indexed variants, biological contexts and probabilities |
| `genome-data` | Noodles VCF ingestion, reference checks and frozen partitions |
| `genome-knowledge` | Evidence, claim transitions and versioned revisions |
| `genome-logic` | Context-sensitive contradictions and implication derivations |
| `genome-proof` | Pinned Lean authority and checked proof artifacts |
| `genome-inference` | Auditable predictor interfaces and simple baselines |
| `genome-agents` | Bounded deterministic research loop and extension interface |
| `genome-evaluation` | Classification metrics, calibration and disagreements |
| `genome-cli` | Demo, mechanism research, study execution and audit |
| `genome-mechanism` | Transcript reconstruction, translation and splice context |
| `genome-atlas` | Frozen primary-chromosome loci and interval coverage |
| `genome-evidence-audit` | Ascertainment, review strength and source/model gaps |
| `genome-research-stats` | Descriptive intervals and paired comparisons |
| `genome-assay` | Reference-mapped reporter effects, units and context |
| `genome-regulatory-data` | Multi-promoter sources, identities, quality and assay roles |
| `genome-regulatory` | Frozen features, ridge fits, guarded tests and controls |

The runtime research agent is deterministic. Live LLM literature retrieval and
foundation-model adapters remain future work. Input types exclude held-out
annotations from predictors, but the process is not a filesystem sandbox for
future shell-capable agents.

Frozen inputs are not silently refreshed during setup. Data-specific cards and
importers describe explicit refresh into new directories; preserve provenance
and revise protocols before treating refreshed data as validation.

See [architecture](docs/architecture.md) and [next scientific goal](research/GOAL.md).
Code is MIT licensed. Datasets, literature and retained reference code keep their
source-specific terms and attribution.
