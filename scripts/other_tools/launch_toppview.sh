#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TOPPVIEW_BIN="$ROOT_DIR/tools/openms/openms-3.5.0/usr/bin/TOPPView"
OPENMS_LIB="$ROOT_DIR/tools/openms/openms-3.5.0/usr/lib"
RUNTIME_LIB="$ROOT_DIR/tools/openms/runtime_root/usr/lib/x86_64-linux-gnu"
QT_PLUGIN_ROOT="$RUNTIME_LIB/qt6/plugins"

if [[ ! -x "$TOPPVIEW_BIN" ]]; then
  echo "TOPPView binary not found: $TOPPVIEW_BIN" >&2
  exit 1
fi

if [[ ! -d "$RUNTIME_LIB" ]]; then
  echo "Missing runtime library directory: $RUNTIME_LIB" >&2
  exit 1
fi

export LD_LIBRARY_PATH="$OPENMS_LIB:$RUNTIME_LIB:${LD_LIBRARY_PATH:-}"
export QT_PLUGIN_PATH="$QT_PLUGIN_ROOT:${QT_PLUGIN_PATH:-}"
export QT_QPA_PLATFORM_PLUGIN_PATH="$QT_PLUGIN_ROOT/platforms"
export QT_QPA_PLATFORM="${QT_QPA_PLATFORM:-xcb}"

if ! "$TOPPVIEW_BIN" "$@"; then
  echo
  echo "TOPPView launch failed." >&2
  echo "Check missing dependencies with:" >&2
  echo "  ldd $TOPPVIEW_BIN | grep 'not found'" >&2
  exit 1
fi
