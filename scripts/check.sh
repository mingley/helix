#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env.sh"
cd "$SCRIPT_DIR/.."
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked --jobs "${HELIX_BUILD_JOBS:-4}" -- -D warnings
cargo test --workspace --locked --jobs "${HELIX_BUILD_JOBS:-4}"
lean proofs/Helix.lean
lean proofs/GenomeMechanisms.lean
python3 -m py_compile scripts/fetch_demo.py
python3 -m py_compile scripts/fetch_hbb_reference.py scripts/fetch_genome_atlas.py
python3 -m py_compile scripts/fetch_hbb_promoter_assay.py
python3 scripts/fetch_hbb_reference.py datasets/hbb-reference --verify
python3 scripts/fetch_genome_atlas.py --self-test
python3 scripts/fetch_genome_atlas.py --verify datasets/genome-atlas
python3 scripts/fetch_hbb_promoter_assay.py --verify datasets/hbb-promoter-assay --reference datasets/hbb-reference
