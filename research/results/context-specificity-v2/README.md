# Context-specificity evidence snapshot

Complete numerical reproduction: **true**.
Full scientific semantic validation: **false**;
status: `unresolved_source_formula_mismatch`. Overall frozen-protocol success
is false. Rust/Python agreement does not resolve the source ESS formula mismatch.
Historical normalization verified: **false**.

This is computational reanalysis of published Element Specificity Score columns,
not raw activity, independent laboratory replication, or a new biological mechanism.
The source's whole-panel normalization includes held-out elements; it is not
prospective or training-only target preprocessing. Published values remain unchanged.
Historical semantic_endpoint_acceptance:true receipts identify the published columns only; their formula acceptance is superseded by the unresolved source-formula audit. They are preserved as historical evidence.

The bound core reports 55,276 sequences,
165,828 measurements, 2,761 groups,
and contexts HepG2, K562, WTC11. Source code commit:
`daea7ed51dd83fe56ec37b0d3ce2f99039f229cb`; working-tree dirty flag: `False`.

- `core/`: schema-4 Rust fits, predictions, features, resampling, ledger, frozen
  protocol, source commitments, reference validation and explicit semantic limits.
- `independent-python.tar.gz`: deterministic archive of the complete independent
  implementation snapshot, source checks, independently fitted outputs and manifests.
- `comparison/`: complete numerical comparison and its implementation snapshot.
- `figures/`: standalone figures, plotted data and exact parent-result commitment.
- `arithmetic/`: post-hoc source arithmetic diagnostic. A fitted weighted plane
  does not identify historical normalization and never changes benchmark targets.
- `endpoint-audit/` and `history/`, when present: independent endpoint audit and
  retained earlier attempts. Superseded results are evidence, not alternate successes.
- `publication.json`: byte sizes and SHA-256 for every other published file,
  plus machine-readable numerical and scientific statuses.

From the repository root, audit the core and verify the complete package:

```bash
source scripts/research-env.sh
cargo run --release --frozen -- audit research/results/context-specificity-v2/core
python3 scripts/package_context_study.py --verify research/results/context-specificity-v2
```

Replay only the saved numerical comparison, using a new replay directory:

```bash
mkdir -p runs/context-published-replay
tar -xzf research/results/context-specificity-v2/independent-python.tar.gz -C runs/context-published-replay
python3 scripts/replicate_context_study.py \
  --compare-rust research/results/context-specificity-v2/core/study.json.gz \
  --independent-run runs/context-published-replay/independent-python \
  --out runs/context-published-replay/comparison
```

Numerical replay returns 0 when agreement holds. To refit both implementations,
regenerate diagnostics and retain the distinct semantic failure, run
`bash scripts/check_context_study.sh runs/<new-prefix>` (replace the placeholder).
The full workflow is expected to return **2** while source semantic validation
remains unresolved, even when its complete numerical comparison matches.

The raw source dataset and its licenses/receipts remain in repository
`datasets/context-panel/`; this bundle commits to them without duplicating them.
The roughly 800 MB hg38 reference cache is separate, under the configured tools
directory; `bash scripts/setup_reference.sh` acquires and verifies it for full replay.
The snapshot is not a standalone replacement for that repository and cache.

Attribute source workbooks to Agarwal et al., *Massively parallel characterization
of transcriptional regulatory elements*, DOI
[10.1038/s41586-024-08430-9](https://doi.org/10.1038/s41586-024-08430-9), and the
preprint-associated deposit [Zenodo 10558183](https://zenodo.org/records/10558183).
The deposit declares [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/);
retain attribution and the source-version distinction recorded in the manifests.
Source material is not relicensed under the repository's code license. Hashes
establish integrity relative to commitments, not empirical truth or source authorship.
