#!/usr/bin/env bash
# Optional full hg38 reference for exact context-panel rederivation.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/env.sh"
REFERENCE_DIR="$HELIX_TOOLS_DIR/reference"
mkdir -p "$REFERENCE_DIR"
REFERENCE_PATH="$REFERENCE_DIR/hg38.2bit"
REFERENCE_PARTIAL="$REFERENCE_PATH.download"
if [ ! -f "$REFERENCE_PATH" ]; then
  REFERENCE_TRANSFER_OK=false
  for REFERENCE_ATTEMPT in 1 2 3; do
    if curl --fail --show-error --location --retry 2 --continue-at - \
      --output "$REFERENCE_PARTIAL" \
      https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.2bit; then
      REFERENCE_TRANSFER_OK=true
      break
    fi
    printf '%s\n' "Reference transfer attempt $REFERENCE_ATTEMPT incomplete; next attempt resumes retained bytes."
  done
  if [ "$REFERENCE_TRANSFER_OK" != true ]; then
    printf '%s\n' 'Reference transfer failed; retained partial file is not accepted.'
    exit 1
  fi
fi
python3 - "$REFERENCE_PATH" "$REFERENCE_PARTIAL" <<'PY'
import hashlib
import os
from pathlib import Path
import sys
final, partial = map(Path, sys.argv[1:])
source = final if final.exists() else partial
md5, sha = hashlib.md5(), hashlib.sha256()
with source.open('rb') as handle:
    for chunk in iter(lambda: handle.read(1024 * 1024), b''):
        md5.update(chunk)
        sha.update(chunk)
if (source.stat().st_size != 835393456
        or md5.hexdigest() != 'dcc3ea27079aa6dc3f9deccd7275e0f8'
        or sha.hexdigest() != '1f67aaa17a77b327738fe750ab37430a85dddf73c7d9189385ad1839259acec0'):
    raise SystemExit('Frozen hg38 reference size/checksum mismatch; file is not accepted')
if source != final:
    os.replace(source, final)
print('Verified frozen hg38 reference:', final)
PY
