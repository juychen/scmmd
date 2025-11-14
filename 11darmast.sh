#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate r_env
# Rscript 11darmast_new.R --input /data1st2/junyi/output/atac1112/subset/region_nt/HIP_HIP_GABA.h5ad \
# --output /data2st1/junyi/output/atac1112/dar/celltype.L2/MAST \

# Rscript 11darmast_new.R --input /data1st2/junyi/output/atac1112/subset/region_nt/PFC_PFC_GABA.h5ad \
# --output /data2st1/junyi/output/atac1112/dar/celltype.L2/MAST \

Rscript 11darmast_new.R --input /data1st2/junyi/output/atac1112/subset/region_nt/HIP_HIP_Glut.h5ad \
--output /data2st1/junyi/output/atac1112/dar/celltype.L2/MAST --region "" \

Rscript 11darmast_new.R --input /data1st2/junyi/output/atac1112/subset/region_nt/PFC_PFC_Glut.h5ad \
--output /data2st1/junyi/output/atac1112/dar/celltype.L2/MAST --region "" \

# Rscript 11darmast_new.R --input /data1st2/junyi/output/atac1112/subset/region_nt/AMY_AMY_Glut.h5ad \
# --output /data2st1/junyi/output/atac1112/dar/celltype.L2/MAST \

# Rscript 11darmast_new.R --input /data1st2/junyi/output/atac1112/subset/region_nt/AMY_AMY_GABA.h5ad \
# --output /data2st1/junyi/output/atac1112/dar/celltype.L2/MAST \

# Rscript 11darmast_new.R --input /data1st2/junyi/output/atac1112/subset/region_nt/AMY_NN.h5ad \
# --output /data2st1/junyi/output/atac1112/dar/celltype.L2/MAST --region AMY

# Rscript 11darmast_new.R --input /data1st2/junyi/output/atac1112/subset/region_nt/PFC_NN.h5ad \
# --output /data2st1/junyi/output/atac1112/dar/celltype.L2/MAST --region PFC

# Rscript 11darmast_new.R --input /data1st2/junyi/output/atac1112/subset/region_nt/HIP_NN.h5ad \
# --output /data2st1/junyi/output/atac1112/dar/celltype.L2/MAST --region HIP
