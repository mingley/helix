# Helix

A Rust research laboratory for typed genomic claims, auditable evidence,
conditional reasoning, and reproducible prediction experiments.

The current [research objective](research/GOAL.md) is to recover a precise HBB
sequence-to-mechanism model, expose its uncertainty, and connect it to measured
regulatory evidence and a genome-wide annotation context. Six research and
engineering subagents contributed independently checked components and source
reviews. This is a partial model of human genome function, not a claim to have
reverse engineered the genome.

## Explore the research workspace

After the setup steps below:

```bash
source scripts/env.sh
cargo run --locked -- research --out runs/hbb-research
cargo run --locked -- audit runs/hbb-research
```

The research command validates the named GRCh38/RefSeq/UniProt sequences,
reconstructs HBB's 628-base transcript and 147-residue protein, enumerates all
4,824 SNVs over its transcript locus, and recomputes consequences for the 54
ClinVar records. It reports coding-only and all-region agreement, evaluates the
computable portions of five hypotheses, audits evidence quality, and places the
locus among 29,949 curated gene loci across the primary reference chromosomes.
It also integrates 561 source-estimated SNV effects from a published HBB
promoter/5-prime-UTR reporter assay, keeping the log2 expression endpoint and
HEL 92.1.7 cell context distinct from clinical classifications.

The output includes a readable research narrative, a standalone interactive
HTML explorer, full JSON annotations and saturation map, source commitments,
an append-only run ledger and scoped Lean proof. Open the generated research.html
locally; no server or remote scripts are needed.

Start with the [biological synthesis](research/hbb-mechanisms.md),
[falsifiable hypotheses](research/hbb-hypotheses.json),
[evidence audit](research/clinvar-audit.md),
[genome context](research/genome-atlas.md), and
[analysis protocol](experiments/hbb-mechanism-protocol.md).
The [reporter dataset card](datasets/hbb-promoter-assay/README.md) and
[assay protocol](experiments/hbb-promoter-assay-protocol.md) document the exact
target mapping, score-set CC0 license, source regression assumptions, retained
equality records and within-study replication limits.
The original clinical heldout has been explored and is now development data.
Its exact paired comparison gives p = 0.375; the earlier accuracy difference
does not establish superiority or independent validation.

The first vertical slice validates 54 frozen public ClinVar HBB variants against
GRCh38, learns a consequence-based classification model from 29 records, predicts
25 held-out records, checks conditional proofs with Lean 4, and saves an auditable
experiment. No API keys, database service, GPU, or LLM is needed.

**This is a working research foundation, not a completed discovery or clinically
validated predictor.** Clinical classifications are not measurements of molecular
function. Logical consistency does not establish biological truth.

## Run

Prerequisites: Linux x86-64, Bash, a C toolchain, curl, and Python 3. The cloud
environment has these prerequisites. Setup installs Rust 1.90.0 and Lean 4.19.0.

```bash
bash scripts/setup.sh
source scripts/env.sh
bash scripts/check.sh
cargo run --locked -- demo --out runs/first
cargo run --locked -- audit runs/first
```

Run commands from the repository root. Tools live outside the checkout:
`/workspace/.tools` in cloud tasks, your local data directory otherwise.
Set `HELIX_TOOLS_DIR` before sourcing `scripts/env.sh` to choose another location.
Setup needs network access; normal demo, tests and audit use the frozen fixture
and installed tools without data retrieval.

Every run requires a new output directory; existing evidence is never overwritten.

| Artifact | Contents |
| --- | --- |
| `report.json` | Dataset commitments, model, metrics, assumptions, proof status, limits |
| `predictions.json` | Rational probabilities for model and prevalence baseline |
| `evidence.json` | Versioned ClinVar accessions, source URLs and hashes, context |
| `ledger.json` | Hash-chained discovery, support, questions, predictions, evaluation, update |
| `proof.lean` | Conditional theorem and concrete probability obligations |
| `artifacts.json` | SHA-256 commitments for the other artifacts |

`audit` rechecks artifacts, revision ordering and the saved proof. Hashes provide
integrity relative to a commitment, not signatures or authentication. Replacing
all commitments can forge a different run. The report records the Git revision
and dirty-state flag when available.

Missing Lean is an error. For explicitly unverified exploration only:

```bash
cargo run --locked -- demo --out runs/unverified --skip-proof
```

That report says `explicitly_skipped` and has no checked proof artifact.
The agent has a default 1,000-operation record budget, hard-capped at 100,000.
Insufficient `--max-steps` fails.

## Scientific contracts

| Layer | Checks | Does not establish |
| --- | --- | --- |
| Rust types | Assembly separation, chromosome bounds, valid probabilities | Empirical truth |
| Ingestion | Coordinates, FASTA alleles, supported SNVs, source metadata and hashes | Indel normalization or liftover |
| Evidence | Matching, non-predicted support in the same context | Direct experimental measurement |
| Replication | Measurements associated with distinct study IDs | Correctness of input study metadata |
| Lean kernel | Entailment from explicit premises and rules; probability bounds | Truth of scientific premises |
| Evaluation | Frozen predictions versus withheld labels and a common baseline | Independent replication or clinical utility |

`Claim<Hypothesis>` becomes `Claim<Supported>` only with matching evidence.
Curated ClinVar assertions cannot become `Claim<Replicated>`. Formal verification
is recorded separately from empirical support; a proof never upgrades evidence.

Coordinates are zero-based inter-residue intervals internally and one-based in
VCF, following the coordinate semantics needed by GA4GH VRS. Display keys such as
`GRCh38:11:5225592:C>A` are **local keys**, not computed GA4GH identifiers.
Complete VRS normalization and identifier interoperability remain future work.

## Workspace

| Crate | Responsibility |
| --- | --- |
| `genome-types` | Assembly-indexed positions and SNVs, contexts, probabilities |
| `genome-data` | Noodles VCF ingestion, reference checks, partitions, provenance |
| `genome-knowledge` | Evidence registry, claim transitions, versioned revisions |
| `genome-logic` | Context-sensitive contradictions and implication derivations |
| `genome-proof` | Pinned Lean subprocess and checked proof artifacts |
| `genome-inference` | Predictor interface, consequence model, prevalence baseline |
| `genome-agents` | Bounded agent interface and deterministic research loop |
| `genome-evaluation` | Accuracy, Brier, log-loss, calibration and disagreements |
| `genome-cli` | Experiment execution and audit |
| `genome-mechanism` | Reference reconstruction, codons, strand-aware SNVs and splice context |
| `genome-atlas` | Frozen primary-chromosome RefSeq loci and merged interval coverage |
| `genome-evidence-audit` | Source ascertainment, review strength, overlap and post-hoc model gaps |
| `genome-research-stats` | Descriptive Wilson intervals, exact paired comparisons and coverage |
| `genome-assay` | Reference-mapped quantitative reporter effects with units and context |

The agent investigates mixed or sparse consequence groups, retrieves training
evidence from the snapshot, proposes independent validation questions, records
assumptions, freezes predictions, and retains unresolved disagreements.
Opposite claims conflict only for the same variant and context. Mixed labels in
a consequence group are model gaps, not automatically contradictory measurements.

The default agent is deterministic. `ResearchAgent` and `Predictor` are extension
interfaces; live LLM literature retrieval and foundation-model adapters are not
implemented. Held-out annotations are excluded from their input types. This
process is not a security sandbox: future shell-capable agents need an
inaccessible holdout store and must not execute arbitrary Lean programs.

## Benchmark limits and data refresh

See [the dataset card](datasets/clinvar-hbb/README.md) for selection, provenance
and split policy. Both models use identical held-out targets. Reports contain
accuracy, balanced accuracy, Brier score, log-loss, five-bin calibration, a
confusion matrix and individual disagreements. Empty/incomplete evaluations fail.

This tiny, selected, imbalanced fixture shares ClinVar curation. Held-out loci may
share submitters and cohorts. Results are descriptive and do not establish
generalization or discovery. Freeze a fresh independent benchmark before real
model selection; repeated tuning on this fixture would make it development data.

To fetch a **new** snapshot without changing the checked-in fixture:

```bash
python3 scripts/fetch_demo.py datasets/clinvar-hbb-new
cargo run --locked -- demo --dataset datasets/clinvar-hbb-new --out runs/new-snapshot
```

This needs network access and can change counts and outcomes. Replacing the frozen
fixture requires provenance review and updated regression expectations; setup
never silently refreshes it.

See [architecture](docs/architecture.md) and [research milestones](docs/roadmap.md).
`bash scripts/check.sh` runs formatting, strict Clippy, Rust tests, Lean theorems,
and the importer syntax check. Integration tests exercise real CLI runs, Lean,
repeatability, failure paths, artifact tampering and immutable output.

Code is MIT licensed. NCBI data retain their original terms and attribution.
Other frozen datasets and literature projections retain the source-specific
terms and attribution documented in their dataset cards and manifests.

For an optional interaction smoke check of a generated explorer (Node.js required):

```bash
node scripts/check_research_viewer.cjs runs/hbb-research/research.html
```

This checks filters, selection and generated DOM behavior without network access.
It is separate from a browser layout or accessibility review.
