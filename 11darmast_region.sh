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
mkdir -p /data2st1/junyi/output/atac1112/dar/region/MASTNG  
Rscript 11darmast_region.R --input /data2st1/junyi/output/atac1112/subset/PFC_20k.h5ad --output /data2st1/junyi/output/atac1112/dar/region/MASTNG --region PFC
Rscript 11darmast_region.R --input /data2st1/junyi/output/atac1112/subset/AMY_20k.h5ad --output /data2st1/junyi/output/atac1112/dar/region/MASTNG --region AMY
Rscript 11darmast_region.R --input /data2st1/junyi/output/atac1112/subset/HIP_20k.h5ad --output /data2st1/junyi/output/atac1112/dar/region/MASTNG --region HIP