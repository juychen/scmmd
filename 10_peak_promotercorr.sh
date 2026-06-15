#!/bin/bash

source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate snapatac2
cd /home/junyichen/code/scmmd/
python peak_promoter_corr_pseudobulk.py \
  --adata-atac /data1st2/junyi/output/atac1112/subset/region_nt/PFC_PFC_GABA.h5ad \
  --promoter-peak /data2st1/junyi/output/atac1112/3REGIONS_peak_l2_promoter_peak_500kb.csv \
  --promoter-col promoter_id \
  --out-prefix /data2st1/junyi/output/atac1112/cCRE/peakpromoter/PFC_PFC_GABA \
  --chunks 40
