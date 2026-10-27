#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/research-env.sh"
cd "$SCRIPT_DIR/.."
mkdir -p "$HELIX_TOOLS_DIR" "$MPLCONFIGDIR" "$PIP_CACHE_DIR" "$XDG_CACHE_HOME/fontconfig"
if [ ! -x "$HELIX_TOOLS_DIR/research-venv/bin/python" ]; then
  python3 -m venv "$HELIX_TOOLS_DIR/research-venv"
fi
"$HELIX_TOOLS_DIR/research-venv/bin/python" -m pip install \
  --disable-pip-version-check --require-hashes -r requirements/research.lock
"$HELIX_TOOLS_DIR/research-venv/bin/python" -c \
  'import sys,numpy,scipy,matplotlib,openpyxl; print("Research Python:",sys.version.split()[0]); print("NumPy",numpy.__version__,"SciPy",scipy.__version__,"Matplotlib",matplotlib.__version__,"openpyxl",openpyxl.__version__)'
printf '%s\n' 'Activate with: source scripts/research-env.sh'
