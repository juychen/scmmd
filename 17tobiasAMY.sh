#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate tobias
ulimit -n 65534
cd /data2st2/junyi/output/atac0627/tobiasbam
for folder in /data2st2/junyi/output/atac0627/tobiasbam/*AMY*; do
  echo $folder
  sample_name=$(basename $folder)
  echo "Processing sample: $sample_name"
  for file in $folder/*; do
    echo $file
    if [[ $file == *.csv ]]; then
      celltype=$(basename $folder)
      echo "Processing celltype: $celltype"
      out_bam="${file%.csv}.bam"
      ctname=$(basename "${file%.csv}")
      out_fp=$folder/corrected/$ctname\_footprints.bw
      mkdir -p $folder/corrected/$ctname
      # if out_bam already exists, skip
      if [[ -f $out_fp ]]; then
        echo "Footprint exist: $out_fp"
        TOBIAS BINDetect --motifs /data2st1/junyi/scenic/mouse/motif/merged_cluster/direct_key_motifs.jaspar \
        --signals $out_fp \
        --genome /data2st1/junyi/ref/GRCm38.p6.genome.fa \
        --peaks /data2st1/junyi/output/atac0627/cCRE/peak.bed \
        --outdir $folder/corrected/$ctname --cores 32 \
        --split 16 --verbosity 4
      fi
    fi
  done
done

# /home/junyichen/subset-bam_linux --bam /data1st2/hannan_25/data/snATAC_process/snATAC_01_bam/MC25A_PFC/outs/possorted_bam.bam \
# --cell-barcodes /data1st2/junyi/output/atac0627/tobiasbam/MC25A_PFC/PFC_L2-3_IT_Glut_MC.csv \
# --out-bam /data1st2/junyi/output/atac0627/tobiasbam/MC25A_PFC/PFC_L2-3_IT_Glut_MC.bam \
# --cores 32