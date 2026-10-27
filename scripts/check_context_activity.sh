#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/research-env.sh"
cd "$SCRIPT_DIR/.."

if [[ "${1:-}" == "--verify-published" ]]; then
  python3 scripts/package_context_activity.py --verify research/results/context-endpoint-recovery-v1
  python3 scripts/package_context_activity.py --verify research/results/context-activity-v3
  cargo run --release --frozen -- audit research/results/context-activity-v3/core
  exit 0
fi

PREFIX="${1:?Usage: bash scripts/check_context_activity.sh runs/NEW_PREFIX | --verify-published}"
for SUFFIX in rust python comparison figures; do
  if [[ -e "${PREFIX}-${SUFFIX}" ]]; then
    printf 'Preserve previous attempt: %s already exists\n' "${PREFIX}-${SUFFIX}" >&2
    exit 1
  fi
done

python3 scripts/build_context_activity.py \
  --source-export research/context-endpoint-recovery/source/jointdata-source-values.json.gz \
  --source-binding research/context-endpoint-recovery/source/schema-binding.json \
  --source-ledger research/context-endpoint-recovery/source/source-type-missing-formula-ledger.jsonl.gz \
  --source-receipt research/context-endpoint-recovery/source/source-export-receipt.json \
  --recovery-summary research/results/context-endpoint-recovery-v1/primary/summary.json \
  --archive-validation research/context-endpoint-recovery/source/archive-validation/verification.json \
  --extraction-receipt research/context-endpoint-recovery/source/member-extraction.json \
  --out datasets/context-activity --verify
cargo run --release --frozen -- context-activity --out "${PREFIX}-rust"
cargo run --release --frozen -- audit "${PREFIX}-rust"
python3 scripts/replicate_context_activity.py --out "${PREFIX}-python"
python3 scripts/replicate_context_activity.py --compare-rust "${PREFIX}-rust/study.json.gz" \
  --independent-run "${PREFIX}-python" --out "${PREFIX}-comparison"
python3 scripts/plot_context_activity.py --source "${PREFIX}-rust" --out "${PREFIX}-figures"
printf '%s\n' 'Published-activity source gates and complete numerical reproduction passed.'
