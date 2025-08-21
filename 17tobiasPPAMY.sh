#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate tobias_cnv
ulimit -n 65534

cd /data2st2/junyi/output/atac0627/tobiasbam
for folder in /data2st2/junyi/output/atac0627/tobiasbam/MW6*AMY*; do
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
      echo "Output BAM: $out_bam"
      
      # if out_bam already exists, skip
      if [[ -f $out_bam ]]; then
        echo "Output BAM file $out_bam already exists. processing in  $folder "
        mkdir -p $folder/corrected
        #cat "test" > $folder/corrected/$sample_name\_$ctname.txt
        out_bw=$folder/corrected/$ctname\_corrected.bw
        # if not exists out_bw, run TOBIAS ATAC correct 

        if [[ ! -f $out_bw ]]; then
          # DRI run echo
          echo "Running TOBIAS ATACorrect for $out_bam"
          TOBIAS ATACorrect --bam $out_bam \
          --genome /data2st1/junyi/ref/GRCm38.p6.genome.fa \
          --peaks /data2st1/junyi/output/atac0627/cCRE/peak.bed \
          --blacklist /data2st1/junyi/ref/mm10-blacklist.v2.bed \
          --outdir $folder/corrected \
          --cores 64
        fi
        # else print exist
        out_fp=$folder/corrected/$ctname\_footprints.bw
        if [[ ! -f $out_fp ]]; then
          echo "Output BW file $out_fp not exists. calculating..."
          echo "generating footprints $out_fp"
          TOBIAS FootprintScores --signal $folder/corrected/$ctname\_corrected.bw \
          --regions /data2st1/junyi/output/atac0627/cCRE/peak.bed \
          --output $out_fp \
          --cores 64
        fi
      fi
      # Run the subset-bam command
      #printf "Running subset-bam for %s with cell barcodes from %s\n" "$sample_name" "$file"
      #/home/junyichen/subset-bam_linux --bam "/data1st2/hannan_25/data/snATAC_process/snATAC_01_bam/$sample_name/outs/possorted_bam.bam" --cell-barcodes "$file" --out-bam "$out_bam" --cores 32
      # Uncomment the line below if you want to run the command for each file
    #   /home/junyichen/subset-bam_linux --bam /data1st2/hannan_25/data/snATAC_process/snATAC_01_bam/$sample_name/outs/possorted_bam.bam \
    #   --cell-barcodes $file \
    #   --out-bam /data1st2/junyi/output/atac0627/tobiasbam/$celltype/${file%.csv}.bam \
    #   --cores 32
    fi
  done
done

# /home/junyichen/subset-bam_linux --bam /data1st2/hannan_25/data/snATAC_process/snATAC_01_bam/MC25A_PFC/outs/possorted_bam.bam \
# --cell-barcodes /data1st2/junyi/output/atac0627/tobiasbam/MC25A_PFC/PFC_L2-3_IT_Glut_MC.csv \
# --out-bam /data1st2/junyi/output/atac0627/tobiasbam/MC25A_PFC/PFC_L2-3_IT_Glut_MC.bam \
# --cores 32