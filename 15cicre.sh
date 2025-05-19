#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate snapatac2

python 15circepipeline.py --input output/atac0416/3REGIONS_peak.h5ads \
--celltype_column region_nt \
--output /data2st1/junyi/output/atac0416/cicre/region_nt \
--region HIP \
--condition MC 

python 15circepipeline.py --input output/atac0416/3REGIONS_peak.h5ads \
--celltype_column region_nt \
--output /data2st1/junyi/output/atac0416/cicre/region_nt \
--region AMY \
--condition MC \

python 15circepipeline.py --input output/atac0416/3REGIONS_peak.h5ads \
--output /data2st1/junyi/output/atac0416/cicre/region_nt \
--celltype_column region_nt \
--region PFC \
--condition MC

python 15circepipeline.py --input output/atac0416/3REGIONS_peak.h5ads \
--celltype_column region_nt \
--output /data2st1/junyi/output/atac0416/cicre/region_nt \
--region HIP \
--condition MW 

python 15circepipeline.py --input output/atac0416/3REGIONS_peak.h5ads \
--celltype_column region_nt \
--output /data2st1/junyi/output/atac0416/cicre/region_nt \
--region AMY \
--condition MW \

python 15circepipeline.py --input output/atac0416/3REGIONS_peak.h5ads \
--output /data2st1/junyi/output/atac0416/cicre/region_nt \
--celltype_column region_nt \
--region PFC \
--condition MW