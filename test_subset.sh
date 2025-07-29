#!/bin/bash
/home/junyichen/subset-bam_linux --bam /data1st2/hannan_25/data/snATAC_process/snATAC_01_bam/MC25A_PFC/outs/possorted_bam.bam \
--cell-barcodes /data1st2/junyi/output/atac0627/tobiasbam/MC25A_PFC/PFC_L2-3_IT_Glut_MC.csv \
--out-bam /data1st2/junyi/output/atac0627/tobiasbam/MC25A_PFC/PFC_L2-3_IT_Glut_MC.bam \
--cores 32