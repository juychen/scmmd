#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate r_env
# Rscript 11darmast.R --input /data1st2/junyi/output/atac0416/subset/region_nt/PFC_PFC_Glut.h5ad \
# --output /data2st1/junyi/output/atac0526/dar/region_nt/mast \
# --celltype_name PFC_PFC_Glut

Rscript 11darmast.R --input /data1st2/junyi/output/atac0416/subset/region_nt/PFC_PFC_GABA.h5ad \
--output /data2st1/junyi/output/atac0526/dar/region_nt/mast \
--celltype_name PFC_PFC_GABA

