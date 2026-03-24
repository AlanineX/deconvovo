#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

TARGETS=(
  "$ROOT_DIR/output/quicklook/01_20260217_wt_TTR_3h/01_20260217_wt_TTR_3h.tic.html"
  "$ROOT_DIR/output/quicklook/03_20260217_wt_TTR_30min/03_20260217_wt_TTR_30min.tic.html"
  "$ROOT_DIR/output/quicklook_deconv/01_20260217_wt_TTR_3h.deconv/01_20260217_wt_TTR_3h.deconv.tic.html"
)

if ! command -v xdg-open >/dev/null 2>&1; then
  echo "xdg-open not found. Open these files manually:" >&2
  printf '%s\n' "${TARGETS[@]}"
  exit 1
fi

for f in "${TARGETS[@]}"; do
  if [[ -f "$f" ]]; then
    xdg-open "$f" >/dev/null 2>&1 &
    sleep 0.3
  else
    echo "Missing quicklook file: $f" >&2
  fi
done
