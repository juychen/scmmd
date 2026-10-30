#!/bin/bash
# =============================================================================
# Region-level subset loom for the CUMS / CSDS single-nucleus data.
#
#   CUMS_M   CUMS  SUS   male    AMY PFC TH HIP
#   CUMS_F   CUMS  SUS   female  AMY PFC TH HIP
#   CSDS_M   CSDS  CSDS  male    AMY PFC TH HIP STR HY MB
#   CON_M    CUMS  CON   male    AMY PFC TH HIP STR HY MB
#   CON_F    CUMS  CON   female  AMY PFC TH HIP
#
# -> AMY.loom PFC.loom TH.loom HIP.loom  (5 groups)
#    STR.loom HY.loom  MB.loom           (CSDS_M + CON_M)
#
# Usage:
#   bash 11subset_loom_meth.sh --dry-run          # only print the cell counts
#   bash 11subset_loom_meth.sh                    # write every region
#   bash 11subset_loom_meth.sh --regions AMY,STR  # a subset of the regions
#   bash 11subset_loom_meth.sh --max_cells_per_group 20000 --overwrite
#   nohup bash 11subset_loom_meth.sh > region_loom.log 2>&1 &
# =============================================================================
set -eo pipefail

source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate scenicplus
ulimit -n 65534

PROJ_DIR="/home/junyichen/code/scmmd"
OUT="/data1st1/junyi/output/meth_region_loom"

mkdir -p "${OUT}"
cd "${PROJ_DIR}"

python 11subsetloom_meth.py \
    --input CUMS=/data2st2/junyi/output/stg1028/CUMS_4VN/CUMS_4VN.h5ad \
    --input CSDS=/data2st2/junyi/output/stg1028/CSDS_4VN/CSDS_4VN.h5ad \
    --output "${OUT}" \
    --layer counts \
    "$@"
