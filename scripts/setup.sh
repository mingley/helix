#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd -- "$SCRIPT_DIR/.." && pwd)"
source "$SCRIPT_DIR/env.sh"
cd "$REPO_DIR"
command -v curl >/dev/null
command -v cc >/dev/null
command -v python3 >/dev/null
mkdir -p "$HELIX_TOOLS_DIR"
SETUP_TEMP="$(mktemp -d)"
trap 'rm -rf -- "$SETUP_TEMP"' EXIT
if [ ! -x "$CARGO_HOME/bin/rustup" ]; then
  curl --fail --show-error --silent --location https://sh.rustup.rs -o "$SETUP_TEMP/rustup.sh"
  sh "$SETUP_TEMP/rustup.sh" -y --no-modify-path --profile minimal --default-toolchain 1.90.0
fi
rustup toolchain install 1.90.0 --profile minimal --component clippy --component rustfmt
if [ ! -x "$ELAN_HOME/bin/elan" ]; then
  curl --fail --show-error --silent --location https://raw.githubusercontent.com/leanprover/elan/v4.1.2/elan-init.sh -o "$SETUP_TEMP/elan.sh"
  sh "$SETUP_TEMP/elan.sh" -y --no-modify-path --default-toolchain leanprover/lean4:v4.19.0
fi
elan toolchain install leanprover/lean4:v4.19.0
rustc --version
lean --version
cargo fetch --locked
cargo build --workspace --locked --jobs "${HELIX_BUILD_JOBS:-4}"
printf '%s\n' 'Setup complete. Run: source scripts/env.sh && bash scripts/check.sh'
