#!/usr/bin/env bash
# Run chromVAR on DMR-filtered peak matrix.
# Usage:
#   bash run_dmr_chromvar.sh [input.h5ad] [output.h5ad]
# Defaults match the project's standard paths.

set -euo pipefail

INPUT="${1:-/data2st1/junyi/output/atac1112/3REGIONS_dmr_filtered.h5ads}"
OUTPUT="${2:-/data2st1/junyi/output/atac1112/chromvardmr/dmr_chromvar.h5ad}"
GENOME="${GENOME:-/data2st1/junyi/ref/GRCm38.p6.genome.fa}"
JASPAR="${JASPAR:-JASPAR2026}"

# DMR matrix is sparse: relax default thresholds
MIN_CELLS_PER_PEAK="${MIN_CELLS_PER_PEAK:-3}"
MIN_GENES_PER_CELL="${MIN_GENES_PER_CELL:-0}"
MIN_TOTAL_COUNTS="${MIN_TOTAL_COUNTS:-200}"

# Whether to downsample to at most N cells (0 = no downsample)
DOWNSAMPLE="${DOWNSAMPLE:-0}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="${SCRIPT_DIR}/dmr_chromvar.py"

# Optional: activate a specific conda env
CONDA_ENV="${CONDA_ENV:-scenicplus}"
if command -v conda >/dev/null 2>&1; then
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate "${CONDA_ENV}"
fi

mkdir -p "$(dirname "${OUTPUT}")"

EXTRA_ARGS=()
if [[ "${DOWNSAMPLE}" -gt 0 ]]; then
  EXTRA_ARGS+=(--downsample-cells "${DOWNSAMPLE}")
else
  EXTRA_ARGS+=(--no-downsample)
fi

echo "========================================================================"
echo "  chromVAR on DMR-filtered peak matrix"
echo "  input  : ${INPUT}"
echo "  output : ${OUTPUT}"
echo "  genome : ${GENOME}"
echo "  jaspar : ${JASPAR}"
echo "  filters: min_cells_per_peak=${MIN_CELLS_PER_PEAK}  min_total_counts=${MIN_TOTAL_COUNTS}"
echo "========================================================================"

python "${SCRIPT}" \
  "${INPUT}" \
  "${OUTPUT}" \
  --genome-file "${GENOME}" \
  --jaspar-release "${JASPAR}" \
  --min-cells-per-peak "${MIN_CELLS_PER_PEAK}" \
  --min-genes-per-cell "${MIN_GENES_PER_CELL}" \
  --min-total-counts "${MIN_TOTAL_COUNTS}" \
  "${EXTRA_ARGS[@]}"
