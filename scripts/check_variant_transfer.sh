#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/research-env.sh"
cd "$SCRIPT_DIR/.."

if [[ "${1:-}" == "--verify-published" ]]; then
  python3 -S scripts/package_variant_transfer.py --verify research/results/variant-transfer-v1
  python3 -S scripts/package_variant_transfer.py --verify research/results/variant-transfer-v2
  exit 0
fi

HELIX_VARIANT_PREFIX="${1:?Usage: bash scripts/check_variant_transfer.sh runs/NEW_PREFIX | --verify-published}"
for HELIX_VARIANT_SUFFIX in observations I S handoff evaluation independent comparison figures; do
  if [[ -e "${HELIX_VARIANT_PREFIX}-${HELIX_VARIANT_SUFFIX}" ]]; then
    printf 'Preserve existing output: %s\n' "${HELIX_VARIANT_PREFIX}-${HELIX_VARIANT_SUFFIX}" >&2
    exit 1
  fi
done
HELIX_VARIANT_BINDING=research/variant-transfer/binding-v2/binding.json
HELIX_VARIANT_BINDING_SHA=a54badc4645d7ac88c31b7722fbbd89aaefa2b995bc0a02befd15123d83d4fba
HELIX_VARIANT_COHORT=research/variant-transfer/binding-v2/cohort.jsonl.gz
HELIX_VARIANT_PROTOCOL=experiments/variant-transfer-v2.json
HELIX_VARIANT_S2="$HELIX_TOOLS_DIR/variant-transfer-v1/S2.xlsx"
HELIX_VARIANT_S4="$HELIX_TOOLS_DIR/variant-transfer-v1/S4.xlsx"

python3 scripts/build_variant_transfer_observations.py --protocol "$HELIX_VARIANT_PROTOCOL" \
  --binding "$HELIX_VARIANT_BINDING" --expected-binding-sha256 "$HELIX_VARIANT_BINDING_SHA" \
  --source "$HELIX_VARIANT_S4" --cohort "$HELIX_VARIANT_COHORT" \
  --out "${HELIX_VARIANT_PREFIX}-observations"
for HELIX_VARIANT_MODEL in I S; do
  python3 scripts/predict_variant_effects.py --model-artifacts research/results/context-activity-v3/core \
    --expected-artifacts-sha256 2ea476bc5413dc13486dec9f3d7e9162a50a0e47f1d83ce29efabc5f5cc5af9d \
    --pairs "${HELIX_VARIANT_PREFIX}-observations/completecasepairs.json" \
    --model "$HELIX_VARIANT_MODEL" --out "${HELIX_VARIANT_PREFIX}-${HELIX_VARIANT_MODEL}"
done
python3 scripts/assemble_variant_transfer.py --observations "${HELIX_VARIANT_PREFIX}-observations" \
  --prediction-I "${HELIX_VARIANT_PREFIX}-I" --prediction-S "${HELIX_VARIANT_PREFIX}-S" \
  --protocol "$HELIX_VARIANT_PROTOCOL" --binding "$HELIX_VARIANT_BINDING" --out "${HELIX_VARIANT_PREFIX}-handoff"
python3 scripts/evaluate_variant_transfer.py --input "${HELIX_VARIANT_PREFIX}-handoff/normalized-input.json" \
  --pins "${HELIX_VARIANT_PREFIX}-handoff/pins.json" --protocol "$HELIX_VARIANT_PROTOCOL" \
  --binding "$HELIX_VARIANT_BINDING" --prediction-I "${HELIX_VARIANT_PREFIX}-I" \
  --prediction-S "${HELIX_VARIANT_PREFIX}-S" --out "${HELIX_VARIANT_PREFIX}-evaluation"
python3 scripts/replicate_variant_transfer.py --run --protocol "$HELIX_VARIANT_PROTOCOL" \
  --binding "$HELIX_VARIANT_BINDING" --expected-binding-sha256 "$HELIX_VARIANT_BINDING_SHA" \
  --cohort "$HELIX_VARIANT_COHORT" --freeze-commit 156c5676c4f79328153edc517a737a2d83c97a48 \
  --exposure-receipt research/variant-transfer/second-freeze-v2.json \
  --s2-workbook "$HELIX_VARIANT_S2" --s4-workbook "$HELIX_VARIANT_S4" --out "${HELIX_VARIANT_PREFIX}-independent"
python3 scripts/replicate_variant_transfer.py --compare-primary "${HELIX_VARIANT_PREFIX}-evaluation" \
  --independent-run "${HELIX_VARIANT_PREFIX}-independent" \
  --primary-observations "${HELIX_VARIANT_PREFIX}-observations" --prediction-I "${HELIX_VARIANT_PREFIX}-I" \
  --prediction-S "${HELIX_VARIANT_PREFIX}-S" --out "${HELIX_VARIANT_PREFIX}-comparison"
python3 scripts/plot_variant_transfer.py --evaluation "${HELIX_VARIANT_PREFIX}-evaluation" \
  --out "${HELIX_VARIANT_PREFIX}-figures"
printf '%s\n' 'Exploratory numerical replay passed. The observed useful-transfer rule is false; retain the negative result.'
