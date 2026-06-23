"""Waters IM-MS analysis pipeline with interactive 2D viewer."""
import os as _os

# Default to a headless matplotlib backend so CLI invocations never pop a window.
# The GUI explicitly switches to its own backend after PySide6 is loaded.
_os.environ.setdefault("MPLBACKEND", "Agg")

__version__ = "0.2.1"
