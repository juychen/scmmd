#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate scenicplus
# cd /data1st2/junyi/output/atac0416/scenic/AMY_AMY_Glut/Snakemake
# # snakemake --cores 32
# cd /data1st2/junyi/output/atac0416/scenic/AMY_AMY_Glut/Snakemake
# snakemake --cores 32

# cd /data1st2/junyi/output/atac0416/scenic/AMY_AMY_Glut/Snakemake
# snakemake --cores 32
cd /data1st2/junyi/output/atac0416/scenic/PFC_PFC_GABA/Snakemake
snakemake --cores 32

cd /data1st2/junyi/output/atac0416/scenic/AMY_AMY_GABA/Snakemake
snakemake --cores 32

cd /data1st2/junyi/output/atac0416/scenic/HIP_HIP_GABA/Snakemake
snakemake --cores 32
