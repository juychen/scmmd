#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate snapatac2

python 15circepipeline.py --input /data1st2/junyi/output/atac1112/subset/region_nt/HIP_HIP_GABA.h5ad --region HIP &\
python 15circepipeline.py --input /data1st2/junyi/output/atac1112/subset/region_nt/HIP_HIP_Glut.h5ad --region HIP ;
python 15circepipeline.py --input /data1st2/junyi/output/atac1112/subset/region_nt/PFC_PFC_GABA.h5ad --region PFC &\
python 15circepipeline.py --input /data1st2/junyi/output/atac1112/subset/region_nt/PFC_PFC_Glut.h5ad --region PFC ;
python 15circepipeline.py --input /data1st2/junyi/output/atac1112/subset/region_nt/AMY_AMY_Glut.h5ad --region AMY &\ 
python 15circepipeline.py --input /data1st2/junyi/output/atac1112/subset/region_nt/AMY_AMY_GABA.h5ad --region AMY ;
python 15circepipeline.py --input /data1st2/junyi/output/atac1112/subset/region_nt/HIP_NN.h5ad --region HIP &\
python 15circepipeline.py --input /data1st2/junyi/output/atac1112/subset/region_nt/PFC_NN.h5ad --region PFC &\
python 15circepipeline.py --input /data1st2/junyi/output/atac1112/subset/region_nt/AMY_NN.h5ad --region AMY ;