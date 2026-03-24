#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="${1:-.}"

mkdir -p \
  "$ROOT_DIR/datasets/incoming/thermo" \
  "$ROOT_DIR/datasets/incoming/waters" \
  "$ROOT_DIR/datasets/incoming/agilent" \
  "$ROOT_DIR/datasets/registered" \
  "$ROOT_DIR/registry" \
  "$ROOT_DIR/standards/mzml" \
  "$ROOT_DIR/standards/im" \
  "$ROOT_DIR/standards/tabular" \
  "$ROOT_DIR/analysis/intact" \
  "$ROOT_DIR/analysis/fragment" \
  "$ROOT_DIR/analysis/hdx" \
  "$ROOT_DIR/analysis/ion-mobility" \
  "$ROOT_DIR/reports"

printf 'Initialized MS workspace scaffold at %s\n' "$(cd "$ROOT_DIR" && pwd)"
