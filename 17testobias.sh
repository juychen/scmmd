#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate snapatac2

# TOBIAS ATACorrect --bam /data1st2/junyi/output/atac0627/tobiasbam/MC25A_PFC/PFC_L2-3_IT_Glut_MC.bam \
# --genome /data2st1/junyi/ref/GRCm38.p6.genome.fa \
# --peaks /data2st1/junyi/output/atac0627/cCRE/peak.bed \
# --blacklist /data2st1/junyi/ref/mm10-blacklist.v2.bed \
# --outdir /data2st1/junyi/output/test \
# --cores 16

# TOBIAS FootprintScores --signal /data2st1/junyi/output/test/PFC_L2-3_IT_Glut_MC_corrected.bw \
# --regions /data2st1/junyi/output/atac0627/cCRE/peak.bed \
# --output /data2st1/junyi/output/test/PFC_L2-3_IT_Glut_MC_footprints.bw --cores 32

TOBIAS BINDetect --motifs /data2st1/junyi/scenic/mouse/motif/merged_cluster/direct_key_motifs.jaspar \
--signals /data2st1/junyi/output/test/PFC_L2-3_IT_Glut_MC_footprints.bw \
--genome /data2st1/junyi/ref/GRCm38.p6.genome.fa \
--peaks /data2st1/junyi/output/atac0627/cCRE/peak.bed \
--outdir /data2st1/junyi/output/test --cond_names PFC_L2-3_IT_Glut_MC --cores 32
