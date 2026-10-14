# -*- mode: python ; coding: utf-8 -*-
"""Freeze only the portable Waters application, including its offline assets."""
import os
from pathlib import Path

import plotly

ROOT = Path(SPECPATH).parent
datas = [
    (str(ROOT / "bin"), "bin"),
    (str(ROOT / "config"), "config"),
    (str(ROOT / "docs"), "docs"),
    (str(ROOT / "README.md"), "."),
    (str(ROOT / "examples" / "production"), "examples/production"),
    (str(Path(plotly.__file__).parent / "package_data"), "plotly/package_data"),
    (os.environ["WATERS_CONSOLE_PROVENANCE"], "build"),
]

a = Analysis(
    [str(ROOT / "build" / "launcher.py")],
    pathex=[str(ROOT), str(ROOT / "build")],
    binaries=[],
    datas=datas,
    hiddenimports=[
        "gui.app", "gui.widgets.waters_viewer", "waters_convert",
        "deconvovo.waters_workflow", "deconvovo.imms_ccs_calibrate",
        "deconvovo.imms_plot", "deconvovo.imms_html",
        "PySide6.QtWebEngineWidgets", "PySide6.QtWebEngineCore",
        "PySide6.QtWebChannel", "plotly.offline",
    ],
    hookspath=[],
    hooksconfig={"matplotlib": {"backends": ["Agg", "QtAgg"]}},
    runtime_hooks=[],
    excludes=[
        "PyQt5", "PyQt6", "PySide2", "tkinter", "pytest", "IPython",
        "torch", "torchvision", "torchaudio", "nvidia", "triton",
        "tensorflow", "onnxruntime", "numba", "llvmlite", "unidec", "pyarrow", "polars",
        "duckdb", "altair", "kaleido", "sqlalchemy", "statsmodels",
        "skimage", "sklearn", "av", "ms_deisotope", "brainpy", "psims",
    ],
    noarchive=False,
)
pyz = PYZ(a.pure)
exe = EXE(
    pyz, a.scripts, [],
    exclude_binaries=True,
    name="WatersConsole",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    disable_windowed_traceback=False,
)
coll = COLLECT(
    exe, a.binaries, a.datas,
    strip=False, upx=False,
    name="WatersConsole",
)
