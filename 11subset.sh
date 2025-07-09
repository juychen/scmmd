#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh

conda activate snapatac2

python 11subseth5.py --input output/atac0627/3REGIONS_peak.h5ads \
--celltype_column region_nt \
--region HIP

python 11subseth5.py --input output/atac0627/3REGIONS_peak.h5ads \
--celltype_column region_nt \
--region AMY

python 11subseth5.py --input output/atac0627/3REGIONS_peak.h5ads \
--celltype_column region_nt \
--region PFC

