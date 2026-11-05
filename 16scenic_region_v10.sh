#!/bin/bash
# =============================================================================
# PySCENIC CTX -> AUCell with the mouse v10nr_clust database
# (region-level looms, sn1029Region)
#
#   GRN is NOT re-run: the v9 adj_<region>.tsv files are reused (GRN does not
#   depend on the motif database).
#
#   input   : /data1st1/junyi/output/sn1029Region/*.loom
#             /data1st1/junyi/output/sn1029Region/scenic_grn/adj_<region>.tsv
#   output  : /data1st1/junyi/output/sn1029Region_v10/scenic_ctx/ctx_<region>.csv
#             /data1st1/junyi/output/sn1029Region_v10/scenic_auc/auc_<region>.tsv
#
# Usage:
#   bash 16scenic_region_v10.sh            # ctx then auc
#   bash 16scenic_region_v10.sh ctx
#   bash 16scenic_region_v10.sh auc
#   nohup bash 16scenic_region_v10.sh > /data1st1/junyi/output/sn1029Region_v10/scenic_v10_run.log 2>&1 &
# =============================================================================
set -uo pipefail

source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate scenicplus
ulimit -n 65534

PROJ_DIR="/home/junyichen/code/scmmd"
INPUT_DIR="/data1st1/junyi/output/sn1029Region"
OUT_DIR="/data1st1/junyi/output/sn1029Region_v10"
ANN="/data2st1/junyi/scenic/mouse/motif/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"
DB="/data2st1/junyi/scenic/mouse/motif/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
PARALLEL=4

STAGE="${1:-all}"

cd "${PROJ_DIR}"
mkdir -p "${OUT_DIR}"

echo "============================================================"
echo " PySCENIC v10nr  |  stage: ${STAGE}"
echo " input : ${INPUT_DIR}"
echo " output: ${OUT_DIR}"
echo " start : $(date '+%F %T')"
echo "============================================================"

if [[ "${STAGE}" == "all" || "${STAGE}" == "ctx" ]]; then
    echo "[$(date '+%F %T')] === CTX (v10nr) ==="
    python 16ctx_region.py \
        --input-dir  "${INPUT_DIR}" \
        --annotations "${ANN}" \
        --database   "${DB}" \
        --output-dir "${OUT_DIR}/scenic_ctx" \
        --parallel   "${PARALLEL}"
    echo "[$(date '+%F %T')] CTX exit=$?"
fi

if [[ "${STAGE}" == "all" || "${STAGE}" == "auc" ]]; then
    echo "[$(date '+%F %T')] === AUCell (v10nr) ==="
    python 16aucell_region.py \
        --input-dir  "${INPUT_DIR}" \
        --ctx-dir    "${OUT_DIR}/scenic_ctx" \
        --output-dir "${OUT_DIR}/scenic_auc" \
        --parallel   "${PARALLEL}"
    echo "[$(date '+%F %T')] AUCell exit=$?"
fi

echo "============================================================"
echo " done  : $(date '+%F %T')"
echo " ctx   : $(ls "${OUT_DIR}"/scenic_ctx/ctx_*.csv 2>/dev/null | wc -l) files"
echo " auc   : $(ls "${OUT_DIR}"/scenic_auc/auc_*.tsv 2>/dev/null | wc -l) files"
echo "============================================================"
