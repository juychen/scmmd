#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate scenicplus

python 15pycistopic.py --input /data1st2/junyi/output/atac0416/subset/region_nt/AMY_AMY_GABA.h5ad \
--output /data1st2/junyi/output/atac0416/cistopic \

# python 15cistopic.py --input /data1st2/junyi/output/atac0416/subset/region_nt/HIP_HIP_GABA.h5ad \
# --output /data1st2/junyi/output/atac0416/cistopic \

# python 15cistopic.py --input /data1st2/junyi/output/atac0416/subset/region_nt/PFC_PFC_GABA.h5ad \
# --output /data1st2/junyi/output/atac0416/cistopic \