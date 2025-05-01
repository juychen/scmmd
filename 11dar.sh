#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh

conda activate snapatac2
python 11darbysubset.py --input /data2st1/junyi/output/atac0416/3REGIONS_peak.h5ads \
--output /data2st1/junyi/output/atac0416/dar \
--region HIP

python 11darbysubset.py --input /data2st1/junyi/output/atac0416/3REGIONS_peak.h5ads \
--output /data2st1/junyi/output/atac0416/dar \
--region AMY

python 11darbysubset.py --input /data2st1/junyi/output/atac0416/3REGIONS_peak.h5ads \
--output /data2st1/junyi/output/atac0416/dar \
--region PFC