# Architecture and trust boundaries

The **17-crate workspace** keeps source binding, numerical reproduction,
conditional formal reasoning and biological validation separate. The current
published-activity benchmark reduces grouped held-out MSE by 8.80% beyond cell
means and passes all 258 independent comparison categories. It predicts the
mean component of three published library-normalized activity scales; it does
not establish absolute expression or a causal regulatory mechanism.

```mermaid
flowchart TD
    W[Official source bytes and exact schema] --> A[Unchanged published RNA/DNA activity]
    W --> E[Published specificity and model-score tables]
    D[Previously reference-validated DNA] --> G[Frozen sequence/locus groups and guards]
    A --> M[Five fresh fixed models]
    G --> M
    M --> N[Independent raw-source, feature, fit and loss replay]
    N --> P[258 numerical comparison categories pass]
    E --> S[Separately frozen arithmetic checks]
    S --> U[Table relationship verified; upstream history separate]
    F[Lean ideal-operator assumptions] --> L[Common activity is not uniquely identifiable from specificity]
    P --> O[Report distinct evidence and limitations]
    U --> O
    L --> O
```

The [source recovery](../research/context-endpoint-recovery.md) obtains official
replicate and aggregate activities and verifies a model-table standardization
relationship. It preserves earlier negative formula checks and their frozen
success rules. The next scientific goal is perturbation prediction in a separate
experimental study, with mechanisms requiring discriminating observations.
See the [activity study](../research/context-activity-study.md) and
[research goal](../research/GOAL.md).

The HBB workflow combines the original clinical demonstration with a reference
mechanism engine, a curated atlas and separately produced reporter evidence:

```mermaid
flowchart TD
    R[GRCh38 genomic sequence and RefSeq transcript/protein] --> M[Validated exon splicing and codon mechanism map]
    U[Reviewed UniProt and source-backed biological questions] --> H[Scoped hypotheses and missing evidence]
    C[Frozen ClinVar assertions] --> A[Annotation comparison and evidence audit]
    M --> A
    G[Curated primary-reference gene annotation] --> X[Genome and locus context]
    Q[MaveDB reporter coefficients and exact target sequence] --> E[Measured log2 reporter endpoint with cell context]
    E --> J[Join exact allele identities while retaining different endpoints]
    A --> J
    J --> O[Research report, explorer and revision ledger]
    X --> O
    H --> O
    M --> P[Small Lean definitions and conditional proofs]
    P --> O
```

Reporter records use a separate quantitative type with units, source accession,
cell line and assay context. They are not converted into clinical labels or used
to promote predicted clinical claims. The assay's target-reference orientation
and the gene's transcription orientation are verified separately.

The research audit checks artifact commitments, ordered revisions, report-to-ledger
bindings, duplicate JSON payloads and the scoped proof. It does not authenticate
source identity or establish biological truth. The hypothesis IDs are fixed by
the finite analysis protocol; alternate files must retain its defined ID set.

The original demo path below remains available and its proof format is preserved.

```mermaid
flowchart TD
    S[Frozen ClinVar and reference FASTA] --> D[Hash and provenance validation]
    D --> T[Assembly-indexed SNVs and context]
    T --> A[Training-only research agent]
    A --> K[Supported claims and revision ledger]
    A --> I[Predictor and prevalence baseline]
    A --> L[Explicit assumptions and implication rules]
    L --> P[Pinned Lean authority]
    I --> F[Freeze predictions]
    F --> E[Withheld-label evaluator]
    P --> R[Auditable experiment artifacts]
    E --> R
    K --> R
```

## Claims and evidence

A proposition identifies a local variant key, biological context and clinical
assertion. A matching curated assertion can support a claim, not establish direct
measurement. Predicted evidence cannot support its own prediction.

Replication promotion requires measured evidence from at least two distinct
nonempty study identities. This checks a data contract; review must still
establish independence and measurement quality. The demo claims no replication.

The evidence registry is content-addressed. Every ledger revision commits to its
payload, event kind, sequence and previous hash. The report commits to inputs,
model counts, predictions and proof.

## Formal boundary

`genome-logic` derives only conclusions reachable from supplied premises through
finite unary implications. It generates propositional symbols internally, so
untrusted atom text cannot inject Lean syntax. Rules carry human-readable
scientific assumptions in the report.

`genome-proof` emits a theorem with premises and rules as hypotheses, plus concrete
natural-number probability bounds. It requires Lean 4.19.0, rejects failed checks,
and terminates proof execution after 30 seconds. The artifact names its authority
and exact source commitment.

Saved proofs are checked against the same finite export grammar before audit
invokes Lean. Imports, evaluation commands, extra axioms and admitted proofs are
rejected rather than treated as user-supplied programs.

The theorem does not establish atom meaning, ClinVar truth, exchangeability,
appropriate decision costs, or relationships in nature. A future LLM can propose
typed statements and rules, not arbitrary Lean programs.

## Evaluation and reproducibility

Fitting accepts training examples. Prediction accepts label-free features.
Evaluation requires exactly one prediction per held-out label. Both models'
predictions enter the ledger before evaluation.

Outputs are deterministic for the same inputs, code revision and authority.
JSON floating-point round trips preserve numeric commitments. The artifact
manifest is written last; incomplete output cannot pass audit.

Hashes are not authentication. Manifests and source URLs need review. A malicious
replacement of all commitments creates another internally consistent experiment.
The process is trusted to follow API boundaries; shell-capable agents require
separate capabilities and an inaccessible holdout store.

## Fixed external transfer and source-gate separation

The variant-transfer path reuses five saved activity fits. Its typed JSON handoff
contains literal paired source effects, exact 200nt alleles, component IDs and
SHA-bound predictions; it accepts no new training or calibration. The observation
builder verifies a committed source/header/cohort binding before parsing selected
S4 cells. Exact feature numerators bridge DNA to fixed fold/ensemble predictions.

V1 retains a strict failed source gate on literal NA. The distinct v2 protocol
declares an exploratory numeric-complete-case population after that exposure,
preserving unknown missingness semantics. Its metadata, preprocessing and numerical
evaluation gates have separate statuses. A useful-transfer result is also separate
from numerical reproduction: all 17 comparison categories pass while the primary
observed improvement is negative and its conditional range crosses zero.

Global grouping is reconstructed by separate Python and C++ implementations on
the complete mapped DNA population before exclusion. The source/cohort seal
preserves original groups through prior-exposure, technical-control and source-ID
filtering. The independent numerical reader separately parses original XLSX DNA
and effects, then rederives features, predictions, losses and resampling with
NumPy. Assembly and source-class preparation remain authenticated handoffs rather
than fully independently rediscovered inputs. Publication manifests retain both
failed and exploratory runs, immutable implementations, ledgers and figures.

`bash scripts/check_variant_transfer.sh runs/NEW_PREFIX` executes the exploratory
numeric replay on the locked cohort; `--verify-published` checks retained bytes
without raw-source access. See the [report](../research/variant-transfer-v2-study.md).

## Current boundaries

The new `helix context-activity` path shares the generic ridge engine but binds a
distinct activity protocol and source. It requires the complete fixed 55,276-ID
cohort, exact official source bytes, metadata-only DNA/group commitments, ordinary
finite numeric cells and jointly bound recovery receipts. Source ESS fields never
enter its response vectors. Schema-5 artifacts retain all fits, observations,
predictions, resamples, feature fractions and source bindings. The audit replays
saved fits and verifies loss and common-component identities; independent Python
separately parses the workbook and refits each model. Its 258-category numerical
comparison is distinct from biological replication.

The generic model's inherited `interpretation_limits` sentence describing primary
group weighting is not the caller's v3 weighting contract: the frozen v3 primary
is equal elements and equal cells, with macro-group summaries as sensitivity
analyses. Actual losses and resampling use the declared equal-element rule.

The `helix context-study` path uses `genome-context-model` for a small generic
matched-context ridge model. The crate enforces engineering contracts; the
frozen scientific protocol separately binds source, endpoint, core, grouping,
guard, features, penalty and interpretation. The CLI verifies the protocol and
input commitments and retains all model fits, predictions, source-feature
dumps, losses and resampling artifacts. The model does not choose a biological
endpoint or establish one from a label's name.

Independent Python reparses the raw workbooks, extracts source230[15:215]
cores, validates every declared hg38 interval/strand, rebuilds group/guard
membership and computes the features without importing Rust. All 55,276
integer-count vectors and 1,824,108 feature fractions match exactly; 829,140
model predictions and every required numerical comparison pass within frozen
tolerances. Independent implementations remain computation checks rather than
independent experimental measurements.

Source ESS semantics have a separate failed/unresolved flag. The observed
columns do not obey the cited zero-sum formula; fitting a weighted plane does
not validate an undocumented normalization. Published values remain unchanged.
The full-library source transformation also prevents a prospective,
training-only target-preprocessing claim. Model training uses only its training
features/labels, but cannot undo the source's panel-level processing.

The complete `scripts/check_context_study.sh` replay produces numerical and
source diagnostics before returning 2 for unresolved scientific validation.
Direct model execution or an integrity audit must not be interpreted as passing
that separate gate. The [context snapshot](../research/results/context-specificity-v2/README.md)
retains the failed source interpretation alongside reproduced computation.

The separate `helix study` path uses `genome-regulatory-data` to validate a
metadata-selected, source-frozen cohort and `genome-regulatory` to execute the
fixed six-model protocol. Model inputs retain sequence and mutation identity;
quality/source fields do not enter the predictor. Training-only scaling,
guarded HBB folds and held-out promoter folds have explicit memberships.
The schema3 artifact includes complete compressed predictions/fits, controls,
context comparisons, input commitments and a revision ledger. Its audit verifies
internal bindings and the scoped method proof; it does not refit every model.

An independent NumPy/SciPy implementation rebuilds features and fits from raw
sources. Continuous predictions, primary losses and controls reproduce, while
26 secondary exact-rank comparisons remain failed. The separately specified
post-hoc finite-grid study is retained alongside that failure. This distinction
is visible in machine-readable flags and the workflow's exit status. Read the
[complete regulatory study](../research/regulatory-study.md) for the scientific
result and limits. The proof covers guard arithmetic, same-input prediction uniqueness
and finite-control resolution, not model coefficients or biological validity.

Mechanism calculations remain limited to the named HBB GRCh38 transcript and
biallelic SNVs. The atlas supplies source-specific genome-wide annotation context;
it does not compute genome-wide variant effects. No full VRS normalization,
indel/liftover engine, validated causal model, general-purpose proof kernel, live
model API, physical experiments or clinical decision workflow is implemented.
