# Helix

A Rust research laboratory for typed genomic claims, auditable evidence,
conditional reasoning, and reproducible prediction experiments.

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
