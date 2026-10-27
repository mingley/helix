# Helix

Helix connects typed genomic representations, source provenance, conditional
Lean reasoning and quantitative experiments. The Rust workspace contains
**17 crates** and builds inspectable partial models of human genome function.

The current [published-activity study](research/context-activity-study.md) tests
whether fixed DNA composition features predict the common component of reporter
activity across HepG2, K562 and WTC11. Shared DNA prediction reduces held-out MSE
by **8.80%** beyond cell means on **55,276 elements in 2,761 sequence/locus groups**.
Independent Rust and Python calculations pass all **258 comparison categories**,
including all 829,140 predictions and exact DNA feature fractions.

An [exact mutation analysis](research/variant-locality-study.md) now identifies
**745,052 of 32,502,288 interior single-base edit opportunities (2.2923%)**
that leave every feature unchanged in this panel. The model must predict zero
for those pairs. It also cannot represent interactions between substitutions
separated by three or more positions. Twenty-one further Lean theorems and
independent exhaustive calculations establish these representation limits;
they do not establish biological neutrality.

The [external variant-effect test](research/variant-transfer-v2-study.md) now
shows **no observed advantage over zero effect** on 2,153 SNVs in 172 components passing the frozen prior-source
exclusion guards from a separate experiment. Macro MSE is 0.2033556 versus 0.2033201;
the conditional improvement interval includes zero. All 17 independent numeric
comparison categories pass. This analysis is explicitly exploratory after a
retained v1 source-token failure; literal NA is excluded operationally with an
unknown missingness reason.

The [source recovery](research/context-endpoint-recovery.md) retrieves the
previously missing official RNA/DNA activity table. It also reproduces all
166,014 model-score cells by standardizing the published ESS table on the exact
model-table subset. Upstream activity-to-ESS processing remains a separate
unresolved question; the earlier frozen source-formula failure is retained.
The source experiment and selected panel are shared with the earlier benchmark.
This is computational reproduction, not an independent laboratory experiment or
a new regulatory mechanism.

[Lean proofs](research/context-identifiability.md) show why centered specificity
alone cannot uniquely identify common activity, including a counterexample with
identical column means and variances. Formal implications, source arithmetic,
predictive evaluation and biological validation remain distinct.

Start with the [research objective](research/GOAL.md),
[activity study](research/context-activity-study.md),
[retained activity evidence](research/results/context-activity-v3/README.md), and
[endpoint recovery](research/context-endpoint-recovery.md).

## Current evidence

| Result | Scope |
| --- | --- |
| Additive DNA/context MSE 0.665083 versus cell means 0.729254 | 8.80% lower held-out MSE on unchanged published log2(RNA/DNA) aggregates |
| Primary MSE gain 0.0641716; conditional 95% range 0.0598440–0.0687169 | 2,000 paired component resamples; no refitting or genome-wide inference |
| All 258 independent activity comparison categories pass | 829,140 model predictions; maximum difference 6.38e−12 |
| S11 ESS → S12 subset column standardization reproduces all 166,014 cells | The older S9 source is byte-identical to final S12; this establishes a file relationship |
| Earlier 36 activity-to-ESS reconstructions all fail | Bounded arithmetic search; further post-hoc cohort diagnostics are separately identified |
| Specificity cannot uniquely determine common activity | 29 checked Lean theorems; ideal operators and explicit moment-preserving witnesses |

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
bash scripts/check_context_activity.sh runs/activity-review
```

Choose a new output prefix. The activity check verifies source and dataset
bindings, runs Rust and independent Python, retains all five models and 2,000
paired resamples, compares complete results, and generates three figures. It
reuses the previously independently verified, unchanged DNA groups and guards.
No new reference download is needed for this activity replay.

The [earlier specificity benchmark](research/context-specificity-study.md)
remains reproducible with
`bash scripts/check_context_study.sh runs/context-review`. That frozen workflow
returns **2** for its original unresolved source-formula gate; later source
recovery does not retroactively change its target or success rule.

Full rederivation obtains a checksum-verified **835 MB hg38 reference**. Cache
it separately if desired:

```bash
bash scripts/setup_reference.sh
```

To inspect the retained evidence without refitting:

```bash
bash scripts/check_context_activity.sh --verify-published
```

The external transfer evidence can also be verified without raw-source access:

```bash
bash scripts/check_variant_transfer.sh --verify-published
```

With the hash-bound Siraj S2/S4 workbooks cached under
`$HELIX_TOOLS_DIR/variant-transfer-v1`, replay the complete exploratory pipeline
with `bash scripts/check_variant_transfer.sh runs/variant-review`. Choose a new
prefix. This reuses the committed cohort, independently reads the original
sources, compares all values/draws and exports standalone figures. Its numerical
pass preserves the negative useful-transfer result. See the
[study and exposure record](research/variant-transfer-v2-study.md).

Run the activity Rust path directly with a new output directory:

```bash
cargo run --release --frozen -- context-activity --out runs/activity-direct-rust
cargo run --release --frozen -- audit runs/activity-direct-rust
```

See the [activity protocol](experiments/context-activity-v3.json),
[independent activity replay](research/independent-context-activity-replication/README.md),
and [source recovery](research/context-endpoint-recovery.md).
The earlier [v2 protocol](experiments/context-specificity-v2.json) and
[independent specificity replay](research/independent-context-replication/README.md)
remain available with their original endpoints.
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
