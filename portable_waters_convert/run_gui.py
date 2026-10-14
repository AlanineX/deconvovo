#!/usr/bin/env python3
"""Launch the Waters GUI with dataset paths relative to this project folder."""
from pathlib import Path
import os


if __name__ == "__main__":
    os.chdir(Path(__file__).resolve().parent)
    from gui.app import main
    main()
