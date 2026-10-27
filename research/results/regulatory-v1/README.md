# Regulatory study v1: published evidence

This is computational reanalysis of one published MPRA study, not independent
laboratory replication or a new biological mechanism. Read the
[research report](../../regulatory-study.md) for findings, methods and limits.

The core was generated from code commit `04ac4bfbe2dfc2f9ee947a5d14968de6cad4ad5f`;
working-tree dirty flag: `False`.

**Primary predictions/losses/controls agree; full frozen numerical reproduction
remains failed, with 26/240 secondary exact-rank failures.** Every failure is in
`comparison/comparison.json`. The separate post-hoc rank supplement reports all
five declared grids; it does not replace the failed endpoint.

- `core/`: complete compressed Rust predictions/fits, scoped proof, source and
  protocol commitments, losses, fold memberships, controls and revision ledger.
- `independent-python.tar.gz`: full independently fitted outputs, features,
  predictions, resampling, software record and exact implementation snapshot.
- `figures/`: five standalone SVG/PNG figures plus exact CSV data and provenance.
- `rank-sensitivity/`: all levels and all metric pairs, including raw failures.
- `history/` (when present): compressed earlier attempts and their original
  manifests/implementations where captured. Superseded artifacts are evidence
  of the investigation, not alternative successful endpoints.
- `publication.json`: SHA-256/size inventory of every published file, including
  histories. It records the partial-validation status explicitly.

From the repository root, verify integrity and replay the saved comparison:

```bash
source scripts/research-env.sh
python3 scripts/package_regulatory_study.py --verify research/results/regulatory-v1
cargo run --release --frozen -- audit research/results/regulatory-v1/core
mkdir -p runs/published-replay
tar -xzf research/results/regulatory-v1/independent-python.tar.gz -C runs/published-replay
python3 scripts/replicate_regulatory_study.py \
  --compare-rust research/results/regulatory-v1/core/study.json.gz \
  --independent-run runs/published-replay/independent-python \
  --out runs/published-replay/comparison
```

The last command returns **2** for the retained comparison failures. Use a new
output directory for another replay. `bash scripts/check_study.sh runs/<new-prefix>`
refits both implementations and regenerates diagnostics and figures.

Source terms remain attached to the source bundles: MaveDB score sets declare
CC0; JASPAR priors are CC BY4.0. Retained author code/documentation are GPL-3.0
reference evidence, and the original OSF table has separately documented producer
terms. Those sources are not relabeled under Helix's MIT code license. Hashes
establish integrity relative to commitments, not authentication or empirical truth.
