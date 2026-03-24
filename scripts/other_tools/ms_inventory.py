#!/usr/bin/env python3
"""Inventory mixed-vendor MS datasets and classify likely workflow types."""

from __future__ import annotations

import argparse
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd


THERMO_PROTEIN_HINTS = ("ttr", "peptide_mapping", "hdx")
WATERS_SMALL_MOLECULE_HINTS = ("adp", "amac", "edda", "water")
WATERS_PROTEIN_HINTS = ("ubq", "myo", "cytc")
AGILENT_PROTEIN_HINTS = ("ubq", "ubiquitin", "ttr", "transthyretin", "cytc", "myo")


def safe_read_text(path: Path, max_chars: int = 200000) -> str:
    try:
        return path.read_text(encoding="utf-8", errors="ignore")[:max_chars]
    except OSError:
        return ""


def parse_waters_header(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    text = safe_read_text(path)
    for line in text.splitlines():
        if not line.startswith("$$ "):
            continue
        payload = line[3:]
        if ":" not in payload:
            continue
        key, value = payload.split(":", 1)
        out[key.strip()] = value.strip()
    return out


def parse_agilent_sample_info(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    try:
        root = ET.parse(path).getroot()
    except (ET.ParseError, OSError):
        return out

    for field in root.findall("./Field"):
        name = (field.findtext("Name") or "").strip()
        value = (field.findtext("Value") or "").strip()
        if name:
            out[name] = value
    return out


def parse_agilent_method(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    text = safe_read_text(path)

    method_name = re.search(r"<MethodName>([^<]+)</MethodName>", text)
    if method_name:
        out["method_name"] = method_name.group(1).strip()

    instrument = re.search(r"<DisplayName>Q-TOF</DisplayName>", text)
    if instrument:
        out["instrument"] = "Q-TOF"

    ion_source = re.search(r"&lt;Name&gt;Ion Source&lt;/Name&gt;.*?&lt;Value&gt;([^<]+)&lt;/Value&gt;", text, re.S)
    if ion_source:
        out["ion_source"] = ion_source.group(1).strip()

    polarity = re.search(r"&lt;Name&gt;Ion Polarity&lt;/Name&gt;.*?&lt;Value&gt;([^<]+)&lt;/Value&gt;", text, re.S)
    if polarity:
        out["polarity"] = polarity.group(1).strip()

    collision = re.search(r"&lt;Name&gt;Collision Energy&lt;/Name&gt;.*?&lt;Value&gt;([^<]+)&lt;/Value&gt;", text, re.S)
    if collision:
        out["collision_energy"] = collision.group(1).strip()

    return out


def classify_thermo_raw(path: Path) -> dict[str, object]:
    lower = path.name.lower()
    workflow = "unknown"
    analyte = "unknown"
    notes: list[str] = []

    if any(h in lower for h in THERMO_PROTEIN_HINTS):
        analyte = "protein_peptide"
        workflow = "fragmented_protein_hdx"
        notes.append("thermo raw names indicate peptide mapping / time-resolved HDX")

    if "peptide_mapping" in lower:
        notes.append("mapping / undeuterated reference run")

    return {
        "vendor": "thermo",
        "format_family": "thermo_raw_file",
        "instrument_hint": "",
        "sample_name": path.stem,
        "method_name": "",
        "analyte_hint": analyte,
        "workflow_hint": workflow,
        "ion_mobility": False,
        "calibrant": False,
        "transform_support": "supported_local_to_mzml",
        "notes": "; ".join(notes),
    }


def classify_waters_raw(path: Path) -> dict[str, object]:
    header = parse_waters_header(path / "_HEADER.TXT")
    acquired_name = header.get("Acquired Name", path.stem)
    instrument = header.get("Instrument", "")
    sample_desc = header.get("Sample Description", "")
    method_name = header.get("MS Method", "")
    lower = f"{path.name} {acquired_name} {sample_desc}".lower()

    analyte = "unknown"
    workflow = "unknown"
    ims_name_text = f"{path.name} {acquired_name}".lower()
    ion_mobility = any(token in ims_name_text for token in ("ims", "mobility", "twims"))
    calibrant = any(token in lower for token in ("cali", "calibration", "nai"))
    notes: list[str] = []

    if any(token in lower for token in WATERS_SMALL_MOLECULE_HINTS):
        analyte = "small_molecule"
        workflow = "intact_small_molecule_screen"
        notes.append("name/header suggest ADP or additive screen")
    elif any(token in lower for token in WATERS_PROTEIN_HINTS):
        analyte = "intact_protein"
        if "den" in lower:
            workflow = "intact_protein_denaturation_or_tuning"
        else:
            workflow = "intact_protein_tuning"
        notes.append("name/header suggest intact protein tune or denaturation run")

    if ion_mobility:
        notes.append("IMS indicated by filename or sample description")
        if workflow == "unknown":
            workflow = "ion_mobility_run"

    if calibrant:
        notes.append("calibrant or calibration marker detected")
        if workflow == "unknown":
            workflow = "ion_mobility_calibrant"

    return {
        "vendor": "waters",
        "format_family": "waters_masslynx_raw_dir",
        "instrument_hint": instrument,
        "sample_name": acquired_name,
        "method_name": method_name,
        "analyte_hint": analyte,
        "workflow_hint": workflow,
        "ion_mobility": ion_mobility,
        "calibrant": calibrant,
        "transform_support": "metadata_only_local",
        "notes": "; ".join(notes),
    }


def classify_agilent_dir(path: Path) -> dict[str, object]:
    sample_info = parse_agilent_sample_info(path / "AcqData" / "sample_info.xml")
    method = parse_agilent_method(path / "AcqData" / "AcqMethod.xml")
    lower = path.name.lower()

    analyte = "unknown"
    workflow = "unknown"
    notes: list[str] = []

    if any(token in lower for token in AGILENT_PROTEIN_HINTS):
        analyte = "intact_protein"
        if any(token in lower for token in ("unfold", "85c", "60c", "30min", "2h", "21h")):
            workflow = "intact_protein_unfolding"
            notes.append("filename indicates thermal or unfolding series")
        elif "fb" in lower:
            workflow = "fragmented_protein_or_front_back_scan"
            notes.append("filename includes fb tag; likely front/back or fragmentation-related acquisition")
        else:
            workflow = "intact_protein_screen"

    return {
        "vendor": "agilent",
        "format_family": "agilent_d_dir",
        "instrument_hint": method.get("instrument", sample_info.get("InstrumentName", "")),
        "sample_name": sample_info.get("Data File", path.stem),
        "method_name": sample_info.get("Method", method.get("method_name", "")),
        "analyte_hint": analyte,
        "workflow_hint": workflow,
        "ion_mobility": False,
        "calibrant": False,
        "transform_support": "metadata_only_local",
        "notes": "; ".join(notes),
    }


def classify_export_csv(path: Path) -> dict[str, object]:
    lower = path.name.lower()
    text = safe_read_text(path, max_chars=4000).lower()
    analyte = "unknown"
    workflow = "tabular_export"
    notes: list[str] = []

    if lower.startswith("mz_abundance_"):
        analyte = "intact_mass_export"
        workflow = "intact_peak_table"
        notes.append("two-column m/z-abundance export")
    elif lower.startswith("fragments_"):
        analyte = "fragment_export"
        workflow = "fragment_annotation_table"
        notes.append("annotated fragment table export")
    elif "x(thomsons)" in text and "y(counts)" in text:
        analyte = "intact_mass_export"
        workflow = "intact_profile_export"
        notes.append("agilent profile export with x(thomsons)/y(counts)")

    return {
        "vendor": "export",
        "format_family": "csv_export",
        "instrument_hint": "",
        "sample_name": path.stem,
        "method_name": "",
        "analyte_hint": analyte,
        "workflow_hint": workflow,
        "ion_mobility": False,
        "calibrant": False,
        "transform_support": "already_tabular",
        "notes": "; ".join(notes),
    }


def iter_artifacts(root: Path) -> list[Path]:
    artifacts: list[Path] = []
    skip_names = {"manifest.csv", "peptide_targets_sequence_start_end.csv"}

    for path in sorted(root.rglob("*")):
        if path.is_symlink():
            continue
        if path.name in skip_names:
            continue
        if path.is_dir() and path.suffix.lower() == ".raw":
            artifacts.append(path)
            continue
        if path.is_dir() and path.suffix.lower() == ".d":
            artifacts.append(path)
            continue
        if path.is_file() and path.suffix.lower() == ".raw":
            artifacts.append(path)
            continue
        if path.is_file() and path.suffix.lower() == ".csv":
            artifacts.append(path)

    return artifacts


def build_rows(roots: list[Path]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    seen: set[Path] = set()
    cwd = Path.cwd().resolve()

    for root in roots:
        for artifact in iter_artifacts(root):
            if artifact in seen:
                continue
            seen.add(artifact)

            if artifact.is_file() and artifact.suffix.lower() == ".raw":
                meta = classify_thermo_raw(artifact)
            elif artifact.is_dir() and artifact.suffix.lower() == ".raw":
                meta = classify_waters_raw(artifact)
            elif artifact.is_dir() and artifact.suffix.lower() == ".d":
                meta = classify_agilent_dir(artifact)
            elif artifact.is_file() and artifact.suffix.lower() == ".csv":
                meta = classify_export_csv(artifact)
            else:
                continue

            try:
                relative_path = str(artifact.resolve().relative_to(cwd))
            except ValueError:
                relative_path = str(artifact.resolve())

            row = {
                "path": str(artifact.resolve()),
                "relative_path": relative_path,
                "root": str(root),
                **meta,
            }
            rows.append(row)

    return rows


def build_summary_markdown(df: pd.DataFrame) -> str:
    lines: list[str] = []
    lines.append("# MS Inventory Summary")
    lines.append("")
    lines.append(f"- Artifacts inventoried: {len(df)}")
    lines.append(f"- Vendors: {', '.join(sorted(df['vendor'].dropna().unique()))}")
    lines.append("")

    for vendor in sorted(df["vendor"].dropna().unique()):
        sub = df[df["vendor"] == vendor].copy()
        lines.append(f"## {vendor.title()}")
        lines.append("")
        lines.append(f"- Artifacts: {len(sub)}")
        families = ", ".join(sorted(sub["format_family"].dropna().unique()))
        lines.append(f"- Formats: {families}")
        workflows = ", ".join(sorted(sub["workflow_hint"].dropna().unique()))
        lines.append(f"- Workflow hints: {workflows}")
        if bool(sub["ion_mobility"].fillna(False).any()):
            lines.append("- Contains ion-mobility candidate runs")
        if bool(sub["calibrant"].fillna(False).any()):
            lines.append("- Contains calibrant candidate runs")
        lines.append("")

        preview = sub[
            [
                "relative_path",
                "format_family",
                "analyte_hint",
                "workflow_hint",
                "ion_mobility",
                "calibrant",
            ]
        ].head(12)
        lines.append(preview.to_csv(index=False))
        lines.append("")

    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description="Inventory mixed-vendor MS datasets")
    parser.add_argument(
        "roots",
        nargs="*",
        default=["data_1_thermo", "data_2_waters", "data_3_agilent"],
        help="Root directories to inventory",
    )
    parser.add_argument("--out-dir", default="output/00_inventory", help="Output directory")
    args = parser.parse_args()

    roots = [Path(root).resolve() for root in args.roots if Path(root).exists()]
    if not roots:
        raise SystemExit("No existing roots provided")

    rows = build_rows(roots)
    df = pd.DataFrame(rows).sort_values(["vendor", "format_family", "relative_path"]).reset_index(drop=True)

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = out_dir / "ms_inventory.csv"
    json_path = out_dir / "ms_inventory.json"
    md_path = out_dir / "ms_inventory_summary.md"

    df.to_csv(csv_path, index=False)
    json_path.write_text(df.to_json(orient="records", indent=2), encoding="utf-8")
    md_path.write_text(build_summary_markdown(df), encoding="utf-8")

    print(f"Wrote: {csv_path}")
    print(f"Wrote: {json_path}")
    print(f"Wrote: {md_path}")
    print(f"Artifacts inventoried: {len(df)}")


if __name__ == "__main__":
    main()
