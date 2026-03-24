#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

echo "[1/5] Organizing data_1 layout"
mkdir -p data_1/archives data_1/extracted data_1/working

for z in data_1/*.zip; do
  [[ -e "$z" ]] || continue
  mv "$z" data_1/archives/
done

for d in data_1/20260217_TTR_TimeP-*; do
  [[ -d "$d" ]] || continue
  mv "$d" data_1/extracted/
done

echo "[2/5] Rebuilding merged working dataset"
rm -rf data_1/working/merged_20260217_TTR_TimeP
python scripts/merge_hdx_dataset.py \
  --input-root data_1/extracted \
  --out-dir data_1/working/merged_20260217_TTR_TimeP \
  --mode symlink

if [[ -f data_1/peptide_targets_sequence_start_end.csv ]]; then
  src_real="$(readlink -f data_1/peptide_targets_sequence_start_end.csv || true)"
  dst_real="$(readlink -f data_1/working/peptide_targets_sequence_start_end.csv || true)"
  if [[ -z "$dst_real" || "$src_real" != "$dst_real" ]]; then
    mv data_1/peptide_targets_sequence_start_end.csv \
      data_1/working/peptide_targets_sequence_start_end.csv
  fi
fi

rm -rf data_1/merged_20260217_TTR_TimeP
ln -sfn working/merged_20260217_TTR_TimeP data_1/merged_20260217_TTR_TimeP
ln -sfn working/peptide_targets_sequence_start_end.csv data_1/peptide_targets_sequence_start_end.csv

echo "[3/5] Organizing output layout"
mkdir -p output

[[ -d output/mzml_merged && ! -L output/mzml_merged ]] && mv output/mzml_merged output/01_mzml
[[ -d output/hdx_raw_reconstructed && ! -L output/hdx_raw_reconstructed ]] && mv output/hdx_raw_reconstructed output/02_hdx_raw
[[ -d output/quicklook && ! -L output/quicklook ]] && mv output/quicklook output/03_quicklook

mkdir -p output/01_mzml output/02_hdx_raw output/03_quicklook
ln -sfn 01_mzml output/mzml_merged
ln -sfn 02_hdx_raw output/hdx_raw_reconstructed
ln -sfn 03_quicklook output/quicklook

echo "[4/5] Writing structure READMEs"
cat > data_1/README.md <<'EOF'
# data_1 Organization

- `archives/`: immutable source zip files
- `extracted/`: extracted source folders
- `working/`: merged/derived working data
- `merged_20260217_TTR_TimeP`: compatibility symlink to `working/merged_20260217_TTR_TimeP`
- `peptide_targets_sequence_start_end.csv`: compatibility symlink to working target CSV
EOF

cat > output/README.md <<'EOF'
# output Organization

- `01_mzml/`: mzML conversions
- `02_hdx_raw/`: reconstructed HDX outputs
- `03_quicklook/`: scan-level QC artifacts
- compatibility symlinks: `mzml_merged`, `hdx_raw_reconstructed`, `quicklook`
- performance summary: `02_hdx_raw/PERFORMANCE.md`
EOF

echo "[5/5] Done"
find data_1 -maxdepth 2 -print | sort
find output -maxdepth 2 -print | sort
