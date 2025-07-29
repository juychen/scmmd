#!/bin/bash
TOBIAS ATACorrect --bam /data1st2/junyi/output/atac0627/tobiasbam/MC25A_PFC/PFC_L2-3_IT_Glut_MC.bam \
--genome /data2st1/junyi/ref/GRCm38.p6.genome.fa \
--peaks /data2st1/junyi/output/atac0627/cCRE/peak.bed \
--blacklist /data2st1/junyi/ref/mm10-blacklist.v2.bed \
--outdir /data2st1/junyi/output/test \
--cores 16