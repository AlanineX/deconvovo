#!/usr/bin/env bash
set -euo pipefail

URL="https://gnps-lcms.ucsd.edu/"

if command -v xdg-open >/dev/null 2>&1; then
  xdg-open "$URL" >/dev/null 2>&1 &
else
  echo "$URL"
fi
