#!/usr/bin/env bash
# Complete the independent work even when a scientific secondary endpoint fails.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/research-env.sh"
cd "$SCRIPT_DIR/.."
if [ "$#" -ne 1 ]; then
  printf '%s\n' 'Usage: bash scripts/check_study.sh runs/<new-prefix>'
  exit 1
fi
STUDY_PREFIX="$1"
for STUDY_SUFFIX in rust python comparison figures rank-sensitivity; do
  if [ -e "$STUDY_PREFIX-$STUDY_SUFFIX" ]; then
    printf '%s\n' "Output exists: $STUDY_PREFIX-$STUDY_SUFFIX"
    exit 1
  fi
done
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cargo run --release --frozen --jobs "${HELIX_BUILD_JOBS:-4}" -- study --out "$STUDY_PREFIX-rust"
cargo run --release --frozen --jobs "${HELIX_BUILD_JOBS:-4}" -- audit "$STUDY_PREFIX-rust"
python3 scripts/replicate_regulatory_study.py --out "$STUDY_PREFIX-python"
COMPARISON_STATUS=0
python3 scripts/replicate_regulatory_study.py \
  --compare-rust "$STUDY_PREFIX-rust/study.json.gz" \
  --independent-run "$STUDY_PREFIX-python" --out "$STUDY_PREFIX-comparison" || COMPARISON_STATUS=$?
if [ "$COMPARISON_STATUS" -ne 0 ] && [ "$COMPARISON_STATUS" -ne 2 ]; then
  exit "$COMPARISON_STATUS"
fi
python3 scripts/check_study_context.py --study "$STUDY_PREFIX-rust" \
  --independent "$STUDY_PREFIX-python" --comparison "$STUDY_PREFIX-comparison"
python3 scripts/plot_regulatory_study.py --run "$STUDY_PREFIX-rust" --out "$STUDY_PREFIX-figures"
python3 scripts/analyze_rank_sensitivity.py --rust "$STUDY_PREFIX-rust/study.json.gz" \
  --independent "$STUDY_PREFIX-python" --out "$STUDY_PREFIX-rank-sensitivity"
printf '%s\n' 'Primary comparison and all secondary failures are retained in the comparison JSON.'
printf '%s\n' 'Post-hoc rank sensitivity is a separately specified diagnostic; see its protocol and script.'
exit "$COMPARISON_STATUS"
