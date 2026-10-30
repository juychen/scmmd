#!/bin/bash
# =============================================================================
# PySCENIC GRN -> CTX -> AUCell for the region-level looms
# (sn1029Region: AMY PFC TH HIP STR HY MB, one loom per brain region)
#
#   input   : /data1st1/junyi/output/sn1029Region/*.loom
#   outputs : scenic_grn/adj_<region>.tsv
#             scenic_ctx/ctx_<region>.csv
#             scenic_auc/auc_<region>.tsv
#
# Each stage skips finished files, so the script can be re-run / resumed.
#
# Usage:
#   bash 16scenic_region.sh                 # run all three stages in order
#   bash 16scenic_region.sh grn             # only GRN
#   bash 16scenic_region.sh ctx             # only CTX   (needs GRN)
#   bash 16scenic_region.sh auc             # only AUCell (needs CTX)
#   nohup bash 16scenic_region.sh > /data1st1/junyi/output/sn1029Region/scenic_run.log 2>&1 &
# =============================================================================
set -uo pipefail

source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate scenicplus
ulimit -n 65534

PROJ_DIR="/home/junyichen/code/scmmd"
INPUT_DIR="/data1st1/junyi/output/sn1029Region"

cd "${PROJ_DIR}"
mkdir -p "${INPUT_DIR}"

STAGE="${1:-all}"

echo "============================================================"
echo " PySCENIC region-level  |  stage: ${STAGE}"
echo " input : ${INPUT_DIR}"
echo " start : $(date '+%F %T')"
echo " looms : $(ls "${INPUT_DIR}"/*.loom 2>/dev/null | wc -l)"
echo "============================================================"

if [[ "${STAGE}" == "all" || "${STAGE}" == "grn" ]]; then
    echo "[$(date '+%F %T')] === GRN ==="
    python 16grn_region.py
    echo "[$(date '+%F %T')] GRN exit=$?"
fi

if [[ "${STAGE}" == "all" || "${STAGE}" == "ctx" ]]; then
    echo "[$(date '+%F %T')] === CTX ==="
    python 16ctx_region.py
    echo "[$(date '+%F %T')] CTX exit=$?"
fi

if [[ "${STAGE}" == "all" || "${STAGE}" == "auc" ]]; then
    echo "[$(date '+%F %T')] === AUCell ==="
    python 16aucell_region.py
    echo "[$(date '+%F %T')] AUCell exit=$?"
fi

echo "============================================================"
echo " done  : $(date '+%F %T')"
echo " grn   : $(ls "${INPUT_DIR}"/scenic_grn/adj_*.tsv 2>/dev/null | wc -l) files"
echo " ctx   : $(ls "${INPUT_DIR}"/scenic_ctx/ctx_*.csv 2>/dev/null | wc -l) files"
echo " auc   : $(ls "${INPUT_DIR}"/scenic_auc/auc_*.tsv 2>/dev/null | wc -l) files"
echo "============================================================"
