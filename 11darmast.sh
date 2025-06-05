#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate r_env
Rscript 11darmast.R --input /data1st2/junyi/output/atac0416/subset/region_nt/AMY_AMY_Glut.h5ad \
--output /data2st1/junyi/output/atac0526/dar/region_nt/mast \
--celltype_name AMY_AMY_Glut

Rscript 11darmast.R --input /data1st2/junyi/output/atac0416/subset/region_nt/AMY_AMY_GABA.h5ad \
--output /data2st1/junyi/output/atac0526/dar/region_nt/mast \
--celltype_name AMY_AMY_GABA

Rscript 11darmast.R --input /data1st2/junyi/output/atac0416/subset/region_nt/HIP_HIP_GABA.h5ad \
--output /data2st1/junyi/output/atac0526/dar/region_nt/mast \
--celltype_name HIP_HIP_GABA

Rscript 11darmast.R --input /data1st2/junyi/output/atac0416/subset/region_nt/HIP_HIP_Glut.h5ad \
--output /data2st1/junyi/output/atac0526/dar/region_nt/mast \
--celltype_name HIP_HIP_Glut
