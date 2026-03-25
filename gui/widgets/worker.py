"""Background worker thread for long-running pipeline tasks."""
from __future__ import annotations

import os
import traceback
from typing import Callable, Any

from PySide6.QtCore import QThread, Signal


# Flag that pipeline code checks to disable multiprocessing inside GUI
_IN_GUI = False


def set_gui_mode():
    """Call at app startup to signal pipeline code to avoid multiprocessing."""
    global _IN_GUI
    _IN_GUI = True
    os.environ["DECONVOVO_GUI"] = "1"


def in_gui() -> bool:
    return _IN_GUI or os.environ.get("DECONVOVO_GUI") == "1"


class Worker(QThread):
    """Run a callable in a background thread with progress/result signals."""

    finished = Signal(object)   # result
    error = Signal(str)         # error message
    progress = Signal(int, int) # (current, total)
    log = Signal(str, str)      # (message, level)

    def __init__(self, fn: Callable, *args, **kwargs):
        super().__init__()
        self._fn = fn
        self._args = args
        self._kwargs = kwargs

    def run(self):
        try:
            import matplotlib
            matplotlib.use("Agg")
            result = self._fn(*self._args, **self._kwargs)
            self.finished.emit(result)
        except Exception as e:
            self.error.emit(f"{type(e).__name__}: {e}\n{traceback.format_exc()}")
