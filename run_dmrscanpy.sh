#!/bin/bash

# Check if the correct number of arguments is provided
# if [ "$#" -ne 2 ]; then
#     echo "Usage: $0 <input_folder> <output_folder>"
#     exit 1
# fi

# Assign input and output folders
source /home/junyichen/anaconda3/etc/profile.d/conda.sh

conda activate snapatac2
# Check if the input folder
python dmr_scanpy.py