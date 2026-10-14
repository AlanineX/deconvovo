#!/usr/bin/env python3
"""Build a self-contained WatersConsole folder on Windows x64 or Linux x64."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_demos(root: Path) -> dict:
    directory = root / "examples" / "production"
    manifest = json.loads((directory / "demo_manifest.json").read_text(encoding="utf-8"))
    count = 0
    for run in manifest["acquisitions"]:
        raw = directory / run["raw_path"]
        actual = {p.relative_to(raw).as_posix() for p in raw.rglob("*") if p.is_file()}
        expected = {entry["path"] for entry in run["files"]}
        if actual != expected:
            raise RuntimeError(f"Incomplete or changed demo acquisition: {run['id']}")
        for entry in run["files"]:
            path = raw / entry["path"]
            if path.stat().st_size != entry["bytes"] or sha256(path) != entry["sha256"]:
                raise RuntimeError(f"Demo differs from its provenance: {path}")
            count += 1
    if len(manifest["acquisitions"]) != 9:
        raise RuntimeError("This release requires all nine complete production RAW runs")
    return {"raw_count": 9, "files_verified": count, "raw_bytes": manifest["bytes"],
            "manifest_sha256": sha256(directory / "demo_manifest.json")}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "build" / "artifacts")
    args = parser.parse_args()
    if platform.machine().casefold() not in {"amd64", "x86_64", "x64"}:
        parser.error("Build with x64 Python: the supplied Waters reader is x64")
    if sys.platform not in {"win32", "linux"}:
        parser.error("This native packaging workflow currently verifies Windows and Linux")
    destination = args.output_dir.expanduser().resolve()
    destination.mkdir(parents=True, exist_ok=True)
    manifest = destination / "build_provenance.json"
    source_files = sorted({
        *ROOT.glob("*.py"), *ROOT.glob("requirements*.txt"),
        *ROOT.glob("deconvovo/**/*.py"), *ROOT.glob("gui/**/*.py"),
        *ROOT.glob("config/*"), *ROOT.glob("bin/*"),
        *ROOT.glob("build/*.py"), *ROOT.glob("build/*.spec"),
        ROOT / "build" / "requirements.txt",
    })
    information = {
        "built_at_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(), "python": sys.version,
        "dependencies": {name: importlib.metadata.version(name) for name in
                         ("PyInstaller", "pyinstaller-hooks-contrib", "PySide6", "plotly",
                          "numpy", "pandas", "scipy", "matplotlib")},
        "demos": verify_demos(ROOT),
        "source_sha256": {p.relative_to(ROOT).as_posix(): sha256(p)
                          for p in source_files if p.is_file()},
    }
    manifest.write_text(json.dumps(information, indent=2) + "\n", encoding="utf-8")
    environment = os.environ.copy()
    environment["WATERS_CONSOLE_PROVENANCE"] = str(manifest)
    subprocess.run(
        [sys.executable, "-m", "PyInstaller", "--noconfirm", "--clean",
         "--distpath", str(destination / "dist"),
         "--workpath", str(destination / "work"),
         str(ROOT / "build" / "waters_console.spec")],
        cwd=ROOT, env=environment, check=True,
    )
    artifact = destination / "dist" / "WatersConsole"
    suffix = ".exe" if os.name == "nt" else ""
    required = [artifact / f"WatersConsole{suffix}"]
    internal = artifact / "_internal"
    required += [internal / "bin" / name for name in
                 ("CDCReader.exe", "MassLynxRaw.dll", "cdt.dll",
                  "Waters_MassLynxSDK_EULA.txt", "UniDec_LICENSE.txt")]
    required += [internal / "plotly" / "package_data" / "plotly.min.js",
                 internal / "config" / "imms_plot_config.json",
                 internal / "build" / manifest.name]
    absent = [str(p) for p in required if not p.is_file()]
    if absent:
        raise RuntimeError("Native bundle is incomplete: " + ", ".join(absent))
    if not any((internal / "PySide6").rglob(f"QtWebEngineProcess{suffix}")):
        raise RuntimeError("Native bundle is missing QtWebEngineProcess")
    if not any((internal / "PySide6").rglob("qtwebengine_resources*.pak")):
        raise RuntimeError("Native bundle is missing the QtWebEngine resource packs")
    for path in source_files:
        relative = path.relative_to(ROOT).as_posix()
        expected = information["source_sha256"].get(relative)
        if expected is not None and sha256(path) != expected:
            raise RuntimeError(f"Source changed while building: {relative}; build again")
        if relative.startswith(("bin/", "config/")) and sha256(internal / relative) != expected:
            raise RuntimeError(f"Frozen asset changed during collection: {relative}")
    if verify_demos(internal) != information["demos"]:
        raise RuntimeError("Frozen production demos do not match the source corpus")
    print(f"Built {artifact}; verified nine complete RAW runs and offline viewer assets.")
    print(f"Run {artifact / ('WatersConsole' + suffix)} --smoke-test OUTPUT_DIRECTORY to validate the native GUI and reader.")


if __name__ == "__main__":
    main()
