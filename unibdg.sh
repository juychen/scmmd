#!/bin/bash

# Check if the correct number of arguments is provided
# if [ "$#" -ne 2 ]; then
#     echo "Usage: $0 <input_folder> <output_folder>"
#     exit 1
# fi

# Assign input and output folders
input_folder="/data2st1/junyi/snmc_bedgraph/peaks"
source /home/junyichen/anaconda3/etc/profile.d/conda.sh

conda activate allcools
# Check if the input folder

cd $input_folder
bedtools unionbedg -i *.bedGraph -header -names $(ls *.bedGraph | sed 's/.bedGraph//g') > /data2st1/junyi/snmc_bedgraph/merged/merge_peak_unified.bedGraph
