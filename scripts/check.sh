#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env.sh"
cd "$SCRIPT_DIR/.."
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked --jobs "${HELIX_BUILD_JOBS:-4}" -- -D warnings
cargo test --workspace --locked --jobs "${HELIX_BUILD_JOBS:-4}"
lean proofs/Helix.lean
python3 -m py_compile scripts/fetch_demo.py
