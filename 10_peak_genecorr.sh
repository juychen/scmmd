#!/bin/bash

# Check if the correct number of arguments is provided
# if [ "$#" -ne 2 ]; then
#     echo "Usage: $0 <input_folder> <output_folder>"
#     exit 1
# fi

# Assign input and output folders
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate snapatac2
cd /home/junyichen/code/scmmd/
python peak_gene_corr_pseudobulk.py \
  --adata-bulk /data2st2/junyi/output/stg1028/CUMS_4VN/CUMS_4VN_bulk.h5ad \
  --adata-atac /data2st1/junyi/output/atac1112/3REGIONS_peak_l2_aggregated.h5ad \
  --promoter-peak /data2st1/junyi/output/atac1112/3REGIONS_peak_l2_promoter_peak_500kb.csv \
  --out-prefix /data2st1/junyi/output/atac1112/cCRE \
  --chunks 40