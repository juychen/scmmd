#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate scenicplus
cd /data2st1/junyi/output/atac1112/scenic/AMY_AMY_Glut/Snakemake
snakemake --cores 32
cd /data2st1/junyi/output/atac1112/scenic/PFC_PFC_Glut/Snakemake
snakemake --cores 32
cd /data2st1/junyi/output/atac1112/scenic/HIP_HIP_Glut/Snakemake
snakemake --cores 32
cd /data2st1/junyi/output/atac1112/scenic/PFC_PFC_GABA/Snakemake
snakemake --cores 32
cd /data2st1/junyi/output/atac1112/scenic/AMY_AMY_GABA/Snakemake
snakemake --cores 32
# cd /data2st1/junyi/output/atac1112/scenic/HIP_HIP_GABA/Snakemake
# snakemake --cores 32
# cd /data2st1/junyi/output/atac1112/scenic/PFC_NN/Snakemake
# snakemake --cores 32
# cd /data2st1/junyi/output/atac1112/scenic/AMY_NN/Snakemake
# snakemake --cores 32
# cd /data2st1/junyi/output/atac1112/scenic/HIP_NN/Snakemake
# snakemake --cores 32
