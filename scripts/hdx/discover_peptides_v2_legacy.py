#!/usr/bin/env python3
"""Discover pepsin-digest peptides from HDX-MS mapping runs.

Uses ms_deisotope for proper isotope pattern scoring:
- PenalizedMSDeconVFitter: MSDeconV × (1 - |G_test|)
- Averagine + brainpy for theoretical isotope distributions
- Charge state enumeration with dependency graph resolution
- Left-search for correct monoisotopic assignment

Then matches deconvoluted envelopes against in-silico pepsin digest.

Usage:
    python scripts/hdx/discover_peptides.py \
        --mzml output/01_hdx_mzml/01_20260217_TTR_peptide_mapping.mzML \
        --known-peptides data_1_thermo/peptide_targets_sequence_start_end.csv \
        --out-dir output/03_hdx_peptide_discovery
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from pyteomics import mass, mzml

PROTON = 1.007276466812

TTR_SEQUENCE = (
    "GSGPTGTGESKCPLMVKVLDAVRGSPAINVAVHVFRKAAD"
    "DTWEPFASGKTSESGELHGLTTEEE"
    "FVEGIYKVEIDT"
    "KSYWKALGISPFHEHA"
    "EVVFTANDSGPRRYTIA"
    "ALLSPYSYSTTAVVTNPKE"
)


def generate_candidates(sequence: str, min_len=5, max_len=40):
    """Generate all sub-peptides with theoretical neutral masses."""
    peptides = {}
    n = len(sequence)
    for start in range(n):
        for end in range(start + min_len, min(start + max_len + 1, n + 1)):
            pep = sequence[start:end]
            if pep not in peptides:
                peptides[pep] = {
                    "sequence": pep, "start": start + 1, "end": end,
                    "neutral_mass": mass.fast_mass(pep, ion_type="M", charge=0),
                }
    return peptides  # {sequence: {start, end, neutral_mass}}


def deconvolute_scan(mz_arr, int_arr, charge_range=(1, 6), ppm=15.0):
    """Deconvolute one MS1 scan using ms_deisotope."""
    from ms_deisotope import deconvolute_peaks
    peaklist = list(zip(mz_arr.astype(float), int_arr.astype(float)))
    if len(peaklist) < 5:
        return []
    result = deconvolute_peaks(
        peaklist,
        charge_range=charge_range,
        error_tolerance=ppm * 1e-6,
        left_search_limit=3,
        truncate_after=0.95,
        iterations=10,
    )
    return list(result.peak_set)


def search_ms1(mzml_path: Path, candidates: dict,
               charge_range=(1, 6), ppm=15.0,
               mass_tol_da=0.02):
    """Deconvolute all MS1 scans, match against candidate peptides."""
    # Build mass lookup: neutral_mass → [peptide_info, ...]
    mass_index = {}
    for pep, info in candidates.items():
        m = round(info["neutral_mass"], 2)
        mass_index.setdefault(m, []).append(info)

    # Sorted masses for binary search
    all_masses = sorted(mass_index.keys())
    all_masses_arr = np.array(all_masses)

    best_hits = {}  # sequence → best hit dict
    n_scans = 0
    n_envelopes_total = 0

    print("  Deconvoluting MS1 scans...")
    for spec in mzml.MzML(str(mzml_path)):
        if spec.get("ms level", 1) != 1:
            continue
        n_scans += 1
        mz_arr = spec.get("m/z array", np.array([]))
        int_arr = spec.get("intensity array", np.array([]))
        if len(mz_arr) < 10:
            continue

        rt = spec.get("scanList", {}).get("scan", [{}])[0].get("scan start time", 0)

        envelopes = deconvolute_scan(mz_arr, int_arr, charge_range, ppm)
        n_envelopes_total += len(envelopes)

        for env in envelopes:
            if env.score < 10:
                continue
            nm = env.neutral_mass

            # Find matching candidates within tolerance
            idx = np.searchsorted(all_masses_arr, nm - mass_tol_da)
            while idx < len(all_masses_arr) and all_masses_arr[idx] <= nm + mass_tol_da:
                for cand in mass_index[all_masses[idx]]:
                    seq = cand["sequence"]
                    mass_err_da = nm - cand["neutral_mass"]
                    mass_err_ppm = mass_err_da / cand["neutral_mass"] * 1e6

                    if abs(mass_err_ppm) > ppm:
                        idx += 1
                        continue

                    if seq not in best_hits or env.score > best_hits[seq]["score"]:
                        # Extract envelope peak info
                        env_peaks = []
                        for ep in env.envelope:
                            env_peaks.append({
                                "mz": round(float(ep.mz), 4),
                                "intensity": round(float(ep.intensity), 1),
                            })

                        best_hits[seq] = {
                            "sequence": seq,
                            "start": cand["start"],
                            "end": cand["end"],
                            "charge": env.charge,
                            "mz_mono": round(float(env.mz), 4),
                            "neutral_mass": round(float(env.neutral_mass), 4),
                            "theo_mass": round(cand["neutral_mass"], 4),
                            "mass_error_ppm": round(mass_err_ppm, 2),
                            "score": round(float(env.score), 2),
                            "n_peaks": len(env.envelope),
                            "scan": n_scans,
                            "rt_min": round(float(rt), 3),
                            "envelope": env_peaks,
                        }
                idx += 1

        if n_scans % 500 == 0:
            print(f"    {n_scans} scans, {n_envelopes_total} envelopes, {len(best_hits)} peptides matched")

    print(f"    {n_scans} scans, {n_envelopes_total} envelopes total")
    print(f"    {len(best_hits)} unique peptides matched")

    if not best_hits:
        return pd.DataFrame()

    df = pd.DataFrame(best_hits.values())
    return df.sort_values(["start", "end"]).reset_index(drop=True)


def compare_known(discovered, known_csv):
    known = pd.read_csv(known_csv)
    known.columns = [c.strip() for c in known.columns]
    known_seqs = set(known["Sequence"].str.upper())
    found_seqs = set(discovered["sequence"].str.upper())

    rows = []
    for _, row in known.iterrows():
        seq = row["Sequence"].upper()
        rows.append({"sequence": seq, "start": row["Start"], "end": row["End"],
                      "in_known": True, "discovered": seq in found_seqs})
    for _, row in discovered.iterrows():
        seq = row["sequence"].upper()
        if seq not in known_seqs:
            rows.append({"sequence": seq, "start": row["start"], "end": row["end"],
                          "in_known": False, "discovered": True,
                          "score": row["score"], "mass_error_ppm": row["mass_error_ppm"]})
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description="Discover HDX-MS peptides (ms_deisotope)")
    parser.add_argument("--mzml", required=True)
    parser.add_argument("--sequence", default=TTR_SEQUENCE)
    parser.add_argument("--known-peptides", default=None)
    parser.add_argument("--out-dir", default="output/peptide_discovery")
    parser.add_argument("--ppm", type=float, default=15.0)
    parser.add_argument("--min-len", type=int, default=5)
    parser.add_argument("--max-len", type=int, default=40)
    parser.add_argument("--charge-min", type=int, default=1)
    parser.add_argument("--charge-max", type=int, default=6)
    parser.add_argument("--min-score", type=float, default=10.0,
                        help="Minimum ms_deisotope score (PenalizedMSDeconV)")
    args = parser.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    seq = args.sequence

    print(f"=== Peptide Discovery ({len(seq)} residues, ms_deisotope) ===")

    # Step 1: Generate candidates
    print("  Generating candidates...")
    candidates = generate_candidates(seq, args.min_len, args.max_len)
    print(f"    {len(candidates)} unique peptides (len {args.min_len}-{args.max_len})")

    # Step 2: Deconvolute + match
    discovered = search_ms1(
        Path(args.mzml), candidates,
        charge_range=(args.charge_min, args.charge_max),
        ppm=args.ppm,
    )
    if discovered.empty:
        print("  No peptides found"); return

    # Filter by score
    discovered = discovered[discovered["score"] >= args.min_score]
    print(f"  After score filter (>={args.min_score}): {len(discovered)} peptides")

    discovered.to_csv(out / "discovered_peptides.csv", index=False, float_format="%.4f")

    # Step 3: Coverage
    coverage = np.zeros(len(seq), dtype=int)
    for _, row in discovered.iterrows():
        coverage[int(row["start"])-1:int(row["end"])] += 1
    cov_pct = np.count_nonzero(coverage) / len(seq) * 100
    print(f"  Coverage: {cov_pct:.1f}%")
    pd.DataFrame({"residue": range(1, len(seq)+1), "aa": list(seq),
                   "coverage": coverage}).to_csv(out / "residue_coverage.csv", index=False)

    # Step 4: Compare known
    if args.known_peptides and Path(args.known_peptides).exists():
        comp = compare_known(discovered, Path(args.known_peptides))
        comp.to_csv(out / "peptide_comparison.csv", index=False)
        n_known = int(comp["in_known"].sum())
        n_found = int(comp[comp["in_known"] & comp["discovered"]].shape[0])
        n_new = int(comp[~comp["in_known"]].shape[0])
        print(f"\n  Fidelity: {n_found}/{n_known} known found ({n_found/max(n_known,1)*100:.0f}%)")
        print(f"  New peptides: {n_new}")

    # Step 5: Score stats
    print(f"\n  Score stats:")
    print(f"    min={discovered['score'].min():.1f}, median={discovered['score'].median():.1f}, "
          f"max={discovered['score'].max():.1f}")
    print(f"    ppm error: mean={discovered['mass_error_ppm'].mean():.2f}, "
          f"std={discovered['mass_error_ppm'].std():.2f}")

    print(f"  Output: {out}")


if __name__ == "__main__":
    main()
