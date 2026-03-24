#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MZMINE_BIN="$ROOT_DIR/tools/mzmine/mzmine-4.9.14/bin/mzmine"
SESSION_HOME="$ROOT_DIR/.tool_home"

if [[ ! -x "$MZMINE_BIN" ]]; then
  echo "MZmine launcher not found or not executable: $MZMINE_BIN" >&2
  exit 1
fi

mkdir -p "$SESSION_HOME" "$SESSION_HOME/.mzmine"

# Keep MZmine logs/config under the project for reproducibility.
export HOME="$SESSION_HOME"
export JAVA_TOOL_OPTIONS="-Duser.home=$SESSION_HOME ${JAVA_TOOL_OPTIONS:-}"

exec "$MZMINE_BIN" "$@"
