#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate snapatac2

python 15export_bw.py


# python 15pycistopic.py --input /data1st2/junyi/output/atac0627/subset/region_nt/AMY_AMY_GABA.h5ad \
# --output /data1st2/junyi/output/atac0627/cistopic \

# python 15pycistopic.py --input /data1st2/junyi/output/atac0627/subset/region_nt/AMY_AMY_Glut.h5ad \
# --output /data1st2/junyi/output/atac0627/cistopic \

# python 15pycistopic.py --input /data1st2/junyi/output/atac0627/subset/region_nt/HIP_HIP_GABA.h5ad \
# --output /data1st2/junyi/output/atac0627/cistopic \

# python 15pycistopic.py --input /data1st2/junyi/output/atac0627/subset/region_nt/HIP_HIP_Glut.h5ad \
# --output /data1st2/junyi/output/atac0627/cistopic \

# python 15pycistopic.py --input /data1st2/junyi/output/atac0627/subset/region_nt/PFC_PFC_GABA.h5ad \
# --output /data1st2/junyi/output/atac0627/cistopic \

# python 15pycistopic.py --input /data1st2/junyi/output/atac0627/subset/region_nt/PFC_PFC_Glut.h5ad \
# --output /data1st2/junyi/output/atac0627/cistopic \