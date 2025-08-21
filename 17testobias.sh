#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate tobias
export OPENBLAS_NUM_THREADS=64
# /home/junyichen/anaconda3/envs/tobias_cnv/lib/python3.12/site-packages/tobias/tools
# set gmm initialization to kmeans ++ if not runnable
ulimit -n 65534
# TOBIAS ATACorrect --bam /data1st2/junyi/output/atac0627/tobiasbam/MC25A_PFC/PFC_L2-3_IT_Glut_MC.bam \
# --genome /data2st1/junyi/ref/GRCm38.p6.genome.fa \
# --peaks /data2st1/junyi/output/atac0627/cCRE/peak.bed \
# --blacklist /data2st1/junyi/ref/mm10-blacklist.v2.bed \
# --outdir /data2st1/junyi/output/test \
# --cores 16

# TOBIAS FootprintScores --signal /data2st1/junyi/output/test/PFC_L2-3_IT_Glut_MC_corrected.bw \
# --regions /data2st1/junyi/output/atac0627/cCRE/peak.bed \
# --output /data2st1/junyi/output/test/PFC_L2-3_IT_Glut_MC_footprints.bw --cores 32
basepath=/data1st2/junyi/output/atac0627/tobiasbam/
TOBIAS BINDetect --motifs /data2st1/junyi/scenic/mouse/motif/merged_cluster/TEST_motifs.jaspar \
--signals $basepath/MC39C_HIP/corrected/HPF_CA1_Glut_MC_footprints.bw  $basepath/MW47A_HIP/corrected/HPF_CA1_Glut_MW_footprints.bw \
--cond_names MC MW \
--genome /data2st1/junyi/ref/GRCm38.p6.genome.fa \
--peaks /data2st1/junyi/output/atac0627/cCRE/peak.bed \
--outdir /data2st1/junyi/output/testtobias --cores 32 \
--split 16 --verbosity 4 \

TOBIAS PlotAggregate --TFBS /data2st1/junyi/output/atac0627/cCRE/peak.bed \
--signals $basepath/MC39C_HIP/corrected/HPF_CA1_Glut_MC_footprints.bw $basepath/MW47A_HIP/corrected/HPF_CA1_Glut_MW_footprints.bw \
--output /data2st1/junyi/output/testtobias/JUND_footprint_HPF_CA1_Glut_MC_footprints.pdf --share_y both --plot_boundaries --signal-on-x