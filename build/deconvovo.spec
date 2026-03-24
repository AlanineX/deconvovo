# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller spec for DeconVoVo Windows build.

Usage (from project root):
    pyinstaller build/deconvovo.spec

This bundles:
- Python runtime + all pip dependencies
- CDCReader.exe + DLLs from UniDec
- deconvovo library + gui
- config/ directory with example CSVs
"""
import sys
from pathlib import Path

block_cipher = None

# Paths
ROOT = Path(SPECPATH).parent
SRC = ROOT / "deconvovo"
GUI = ROOT / "gui"
CONFIG = ROOT / "config"

# Find CDCReader.exe and DLLs from UniDec
datas = []
try:
    import unidec
    unidec_dir = Path(unidec.__file__).parent
    cdcreader = unidec_dir / "bin" / "CDCReader.exe"
    if cdcreader.exists():
        # Bundle entire bin/ directory (CDCReader + DLLs)
        datas.append((str(unidec_dir / "bin"), "unidec/bin"))
    # Also bundle Waters importer DLLs if present
    waters_dir = unidec_dir / "bin" / "UniDecImporter" / "Waters"
    if waters_dir.exists():
        datas.append((str(waters_dir), "unidec/bin/UniDecImporter/Waters"))
except ImportError:
    print("WARNING: UniDec not found — CDCReader.exe will not be bundled")

# Bundle config/ with example CSVs
if CONFIG.exists():
    datas.append((str(CONFIG), "config"))

a = Analysis(
    [str(ROOT / "gui" / "app.py")],
    pathex=[str(ROOT)],
    binaries=[],
    datas=datas,
    hiddenimports=[
        "deconvovo",
        "deconvovo.io",
        "deconvovo.vendors",
        "deconvovo.vendors.waters",
        "deconvovo.smooth",
        "deconvovo.parallel",
        "deconvovo.taskqueue",
        "deconvovo.ccs_peak_pick",
        "deconvovo.ccs_fit",
        "deconvovo.ccs_convert",
        "deconvovo.ccs_plot",
        "deconvovo.imms_ccs_calibrate",
        "deconvovo.imms_html",
        "deconvovo.imms_plot",
        "deconvovo.imms_deconv",
        "deconvovo.imms_summary",
        "deconvovo.imms_convert",
        "deconvovo.waters_convert",
        "scipy.signal",
        "scipy.ndimage",
        "scipy.stats",
        "plotly",
        "unidec",
        "unidec.engine",
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=["tkinter", "test", "unittest"],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="DeconVoVo",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,  # no terminal window
    disable_windowed_traceback=False,
    argv_emulation=False,
    icon=None,  # TODO: add icon
)

coll = COLLECT(
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name="DeconVoVo",
)
