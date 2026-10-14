"""Frozen entry point; the validation switch is for native build checks only."""
from __future__ import annotations

import os
from pathlib import Path
import sys


def main() -> None:
    root = Path(sys._MEIPASS) if getattr(sys, "frozen", False) else Path(__file__).resolve().parents[1]
    os.chdir(root)
    if not getattr(sys, "frozen", False):
        sys.path.insert(0, str(root))
    if len(sys.argv) == 3 and sys.argv[1] == "--smoke-test":
        from native_smoke import run
        raise SystemExit(run(Path(sys.argv[2]).expanduser().resolve(), root))
    from gui.app import main as gui_main
    gui_main()


if __name__ == "__main__":
    main()
