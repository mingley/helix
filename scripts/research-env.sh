#!/usr/bin/env bash
HELIX_RESEARCH_SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$HELIX_RESEARCH_SCRIPT_DIR/env.sh"
export UV_CACHE_DIR="$HELIX_TOOLS_DIR/uv-cache"
export PIP_CACHE_DIR="$HELIX_TOOLS_DIR/pip-cache"
export MPLCONFIGDIR="$HELIX_TOOLS_DIR/matplotlib-cache"
export XDG_CACHE_HOME="$HELIX_TOOLS_DIR/research-xdg-cache"
# Imported archived implementations must not gain unmanifested cache files.
export PYTHONDONTWRITEBYTECODE=1
if [ -x "$HELIX_TOOLS_DIR/research-venv/bin/python" ]; then
  export PATH="$HELIX_TOOLS_DIR/research-venv/bin:$PATH"
fi
