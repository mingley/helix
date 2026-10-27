#!/usr/bin/env bash
# Complete numerical work and retain the distinct unresolved semantic failure.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/research-env.sh"
cd "$SCRIPT_DIR/.."
if [ "$#" -ne 1 ]; then
  printf '%s\n' 'Usage: bash scripts/check_context_study.sh runs/<new-prefix>'
  exit 1
fi
CONTEXT_PREFIX="$1"
for CONTEXT_SUFFIX in rust python comparison figures arithmetic endpoint-audit; do
  if [ -e "$CONTEXT_PREFIX-$CONTEXT_SUFFIX" ]; then
    printf '%s\n' "Output exists: $CONTEXT_PREFIX-$CONTEXT_SUFFIX"
    exit 1
  fi
done
if [ -e "$CONTEXT_PREFIX-groups.json" ]; then
  printf '%s\n' "Output exists: $CONTEXT_PREFIX-groups.json"
  exit 1
fi
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
bash scripts/setup_reference.sh
python3 scripts/fetch_context_panel.py verify datasets/context-panel \
  --reference "$HELIX_TOOLS_DIR/reference/hg38.2bit"
python3 scripts/group_context_panel.py --panel datasets/context-panel/panel.json \
  --out "$CONTEXT_PREFIX-groups.json"
cargo run --release --frozen --jobs "${HELIX_BUILD_JOBS:-4}" -- context-study \
  --groups "$CONTEXT_PREFIX-groups.json" --out "$CONTEXT_PREFIX-rust"
cargo run --release --frozen --jobs "${HELIX_BUILD_JOBS:-4}" -- audit "$CONTEXT_PREFIX-rust"
python3 scripts/replicate_context_study.py --dataset datasets/context-panel \
  --groups "$CONTEXT_PREFIX-groups.json" --reference "$HELIX_TOOLS_DIR/reference/hg38.2bit" \
  --out "$CONTEXT_PREFIX-python"
CONTEXT_COMPARE_STATUS=0
python3 scripts/replicate_context_study.py --compare-rust "$CONTEXT_PREFIX-rust/study.json.gz" \
  --independent-run "$CONTEXT_PREFIX-python" --out "$CONTEXT_PREFIX-comparison" || CONTEXT_COMPARE_STATUS=$?
if [ "$CONTEXT_COMPARE_STATUS" -ne 0 ] && [ "$CONTEXT_COMPARE_STATUS" -ne 2 ]; then
  exit "$CONTEXT_COMPARE_STATUS"
fi
python3 scripts/replicate_context_study.py --audit-endpoint \
  --independent-run "$CONTEXT_PREFIX-python" --out "$CONTEXT_PREFIX-endpoint-audit"
python3 scripts/audit_context_source_arithmetic.py --out "$CONTEXT_PREFIX-arithmetic"
python3 scripts/plot_context_study.py --run "$CONTEXT_PREFIX-rust" \
  --arithmetic "$CONTEXT_PREFIX-arithmetic" --out "$CONTEXT_PREFIX-figures"
python3 - "$CONTEXT_PREFIX" <<'PY'
import json
from pathlib import Path
import sys
prefix = sys.argv[1]
comparison = json.loads(Path(prefix + '-comparison/comparison.json').read_bytes())
summary = json.loads(Path(prefix + '-rust/summary.json').read_bytes())
print('Full independent numerical reproduction:', comparison['matched'])
semantic = summary['source_semantics']['full_scientific_semantic_validation']
print('Full scientific source-formula validation:', semantic)
print('Published source targets remain unchanged. Numeric agreement does not resolve their formula.')
raise SystemExit(0 if comparison['matched'] and semantic else 2)
PY
