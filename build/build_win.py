#!/usr/bin/env python3
"""Build the verified WatersConsole folder using the pinned portable builder."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent.parent


def main() -> None:
    subprocess.run(
        [sys.executable, str(ROOT / "portable_waters_convert/build/build_native.py"),
         "--output-dir", str(ROOT / "waters-build")],
        cwd=ROOT, check=True,
    )


if __name__ == "__main__":
    main()
