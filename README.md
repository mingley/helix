# Helix

Helix connects typed genomic representations, source provenance, conditional
Lean reasoning and quantitative experiments. The Rust workspace contains
**17 crates** and builds inspectable partial models of human genome function.

The current [context-score benchmark](research/context-specificity-study.md)
reduces held-out MSE by **13.434%** on fixed published source columns. Independent
implementations reproduce every required numerical comparison, including exact
features and all 829,140 model predictions.

**Scientific source-formula validation remains unresolved.** The supplied
columns fail the published ESS zero-sum formula by far more than rounding.
Numerical reproduction does not resolve their preprocessing history or certify
the stated biological endpoint. The workflow retains this failure and returns
exit status 2. No target projection, repair or refit is used to force agreement.

The next goal is to recover authoritative, version-linked normalization/raw
activity or independently compatible measurements. No new biological mechanism
or global novelty is established.

Start with the [research objective](research/GOAL.md),
[context study](research/context-specificity-study.md),
[retained snapshot](research/results/context-specificity-v2/README.md), and
[source-semantic review](research/context-specificity-semantic-review.md).

## Current evidence

| Result | Scope |
| --- | --- |
| Context-specific MSE 0.86538919 versus additive 0.99969020 | 13.434% lower MSE on 55,276 selected, reference-valid elements in three known contexts |
| All 243 independent comparison categories pass | 1,824,108 feature fractions match exactly; 829,140 model predictions reproduce within frozen tolerance |
| Every admitted score triple violates the stated zero-sum formula | Maximum absolute sum 1.0987267; source-normalization provenance remains unresolved |
| Source columns instead occupy a fitted weighted rank-two plane | Compatible with extra column scaling; no proof of historical preprocessing or new biology |
| Raw-activity context-prediction-v1 stopped before fitting | Available files contained ESS rather than the required activity means; v2 is a distinct protocol |

The earlier [regulatory study](research/regulatory-study.md) remains negative:
none of five fixed nonconstant models beats the mean across nine held-out
promoters. Its 26/240 exact-rank failures remain failed; separate post-hoc grid
agreement does not repair that endpoint. HBB deletion/QC forensics remain
source accounting, not a biological noise floor. Known context/motif prediction
has substantial prior art.

## Set up and reproduce the study

Validated environment: Linux x86-64, Rust 1.90.0, Lean 4.19.0 and CPython 3.12.
Bash, curl and a C toolchain are prerequisites. No API keys, database, GPU or live
LLM is required. Run from the repository root:

```bash
bash scripts/setup.sh
bash scripts/setup_research.sh
source scripts/research-env.sh
bash scripts/check.sh
bash scripts/check_context_study.sh runs/context-review
```

Choose a new output prefix. The context check rebuilds groups and guards,
replays reference validation, runs Rust and independent Python, retains all five
models and 2,000 paired resamples, compares artifacts, and generates figures
and source diagnostics. **Exit status 2 is expected because scientific
source-formula validation remains unresolved, despite numerical reproduction
passing.** It must not be presented as a complete scientific pass.

Full rederivation obtains a checksum-verified **835 MB hg38 reference**. Cache
it separately if desired:

```bash
bash scripts/setup_reference.sh
```

Snapshot inspection and direct fitting of the frozen panel do not require that
download. To create a new grouping receipt and run the Rust path directly:

```bash
python3 scripts/group_context_panel.py --panel datasets/context-panel/panel.json \
  --out runs/context-direct-groups.json
cargo run --release --frozen -- context-study \
  --groups runs/context-direct-groups.json --out runs/context-direct-rust
cargo run --release --frozen -- audit runs/context-direct-rust
```

See the [independent replay](research/independent-context-replication/README.md),
[v2 protocol](experiments/context-specificity-v2.json), and
[source interpretation](research/context-literature.md).
The [earlier study snapshot](research/results/regulatory-v1/README.md) remains
available. Reproduce it with `bash scripts/check_study.sh runs/study-review`;
that workflow returns 2 for its retained secondary-rank failures.

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
| `genome-context-model` | Matched-context sequence features, ridge models and paired losses |

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
