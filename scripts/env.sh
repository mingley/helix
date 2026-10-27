#!/usr/bin/env bash
# Source from any directory. Keep tools outside the checkout in cloud tasks.
if [ -z "${HELIX_TOOLS_DIR:-}" ]; then
  if [ -d /workspace ]; then
    export HELIX_TOOLS_DIR=/workspace/.tools
  else
    export HELIX_TOOLS_DIR="${XDG_DATA_HOME:-$HOME/.local/share}/helix"
  fi
fi
export CARGO_HOME="$HELIX_TOOLS_DIR/cargo"
export RUSTUP_HOME="$HELIX_TOOLS_DIR/rustup"
export ELAN_HOME="$HELIX_TOOLS_DIR/elan"
export PATH="$CARGO_HOME/bin:$ELAN_HOME/bin:$PATH"
