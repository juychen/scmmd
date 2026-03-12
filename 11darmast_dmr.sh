#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate r_env
# Rscript 11darmast_new.R --input /data1st2/junyi/output/atac1112/subset/region_nt/HIP_HIP_GABA.h5ad \
# --output /data2st1/junyi/output/atac1112/dar/region/MAST \

# Rscript 11darmast_new.R --input /data1st2/junyi/output/atac1112/subset/region_nt/PFC_PFC_GABA.h5ad \
# --output /data2st1/junyi/output/atac1112/dar/region/MAST \

# Rscript 11darmast_new.R --input /data1st2/junyi/output/atac1112/subset/region_nt/AMY_AMY_Glut.h5ad \
# --output /data2st1/junyi/output/atac1112/dar/region/MAST \

# Rscript 11darmast_new.R --input /data1st2/junyi/output/atac1112/subset/region_nt/AMY_AMY_GABA.h5ad \
# --output /data2st1/junyi/output/atac1112/dar/region/MAST \
mkdir -p /data2st1/junyi/output/atac1112/dar/region/MASTdmr  
Rscript 11darmast_new_dmr.R --input /data2st1/junyi/output/atac1112/subset/dmr_region_nt/PFC_5hmC.h5ad --output /data2st1/junyi/output/atac1112/dar/region/MASTdmr --region "PFC" &\