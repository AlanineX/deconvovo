#!/usr/bin/env python3
"""Build DeconVoVo Windows executable.

Run on Windows:
    python build/build_win.py

Requires: pip install pyinstaller
Produces: dist/DeconVoVo/ folder ready to zip and distribute.
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
SPEC = ROOT / "build" / "deconvovo.spec"


def main():
    print("Building DeconVoVo Windows executable...")
    print(f"  Spec: {SPEC}")
    print(f"  Root: {ROOT}")

    cmd = [
        sys.executable, "-m", "PyInstaller",
        "--noconfirm",
        "--clean",
        str(SPEC),
    ]

    print(f"  Command: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=str(ROOT))

    if result.returncode == 0:
        dist = ROOT / "dist" / "DeconVoVo"
        print(f"\nBuild successful!")
        print(f"  Output: {dist}")
        print(f"  To distribute: zip the {dist} folder")
    else:
        print(f"\nBuild failed with exit code {result.returncode}")
        sys.exit(1)


if __name__ == "__main__":
    main()
