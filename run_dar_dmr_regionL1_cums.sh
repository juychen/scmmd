#!/bin/bash
# =============================================================================
# DAR on DMR regions -- REGION x (celltype.L1, Neuron split into GABA/Glut)
#                        CUMS vs CONTROL, MAST + NG covariate
#
#   grouping  : region_l1  (18 groups)
#                 AMY_GABA, AMY_Glut, AMY_Astro-Epen, AMY_Immune,
#                 AMY_OPC-Oligo, AMY_Vascular
#                 HIP_... (same six)
#                 PFC_... (same six)
#               neurons are split with celltype.L1_ct (GABA / Glut);
#               non-neurons keep their celltype.L1 label.
#   contrast  : CONTROL (MW)  ->  CUMS (MC)      logFC positive = higher in CUMS
#   input     : 3REGIONS_dmr_filtered.h5ads      (74,477 DMRs x 152,117 cells)
#   model     : zlm(~ compare_group + ngeneson)
#               - normalised with obs$total_counts  (FULL peak matrix depth)
#               - ngeneson taken from obs            (FULL peak matrix, z-scored)
#
# Usage:
#   bash run_dar_dmr_regionL1_cums.sh
#   nohup bash run_dar_dmr_regionL1_cums.sh > regionL1_cums.nohup.log 2>&1 &
# =============================================================================
set -eo pipefail

source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate r_env

## ---------------------------------------------------------------- settings --
PROJ_DIR="/home/junyichen/code/scmmd"
INPUT="/data2st1/junyi/output/atac1112/3REGIONS_dmr_filtered.h5ads"
OUTDIR="/data2st1/junyi/output/atac1112/dar/dmr_region_l1_cums"
QC_CSV="/data2st1/junyi/output/atac0627/frac_qc.csv"

R_SCRIPT="${PROJ_DIR}/11darmast_dmr_cums.R"

# composite grouping: Region x celltype.L1, neurons split into GABA / Glut
GROUP_BY="region_l1"

# contrast: group1 = CONTROL, group2 = CUMS   (MW = CON, MC = SUS in this dataset)
COMPARE_BY="Condition"
GROUP1="MW"     # CONTROL
GROUP2="MC"     # CUMS / susceptible

MAX_CELLS=10000

## ------------------------------------------------------------------- checks --
[[ -f "${R_SCRIPT}" ]] || { echo "[ERROR] R script not found: ${R_SCRIPT}"; exit 1; }
[[ -f "${INPUT}"    ]] || { echo "[ERROR] input h5ad not found: ${INPUT}";   exit 1; }
if [[ ! -f "${QC_CSV}" ]]; then
    echo "[WARN] QC csv not found (${QC_CSV}); batch/frac columns will be skipped"
fi

mkdir -p "${OUTDIR}"
cd "${PROJ_DIR}"

echo "============================================================"
echo " job      : DAR on DMR / Region x L1 (GABA,Glut) / ${GROUP1}(CONTROL) vs ${GROUP2}(CUMS)"
echo " input    : ${INPUT}"
echo " output   : ${OUTDIR}"
echo " start    : $(date '+%F %T')"
echo "============================================================"

## --------------------------------------------------------------------- run --
Rscript "${R_SCRIPT}" \
    --input                "${INPUT}" \
    --output               "${OUTDIR}" \
    --group_by             "${GROUP_BY}" \
    --compare_by           "${COMPARE_BY}" \
    --group1               "${GROUP1}" \
    --group2               "${GROUP2}" \
    --batch_by             "date" \
    --use_batch            "no" \
    --exclude              "" \
    --prefix               "regionL1_cums" \
    --qc_csv               "${QC_CSV}" \
    --max_cells            "${MAX_CELLS}" \
    --subsample_threshold  25000 \
    --min_cells_peak       10 \
    --freq_expressed       0.1 \
    --scale_factor         100000 \
    2>&1 | tee "${OUTDIR}/run_dar_dmr_regionL1_cums.log"

echo "============================================================"
echo " done     : $(date '+%F %T')"
echo " results  : ${OUTDIR}/MASTNG_regionL1_cums_all_comparisons.csv"
echo "============================================================"
