#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh

conda activate snapatac2

python 11darbysubset.py --input output/atac0416/3REGIONS_peak.h5ads \
--output /data2st1/junyi/output/atac0526/dar \
--celltype_column region_nt \
--method memento-binary \
--region PFC

python 11darbysubset.py --input output/atac0416/3REGIONS_peak.h5ads \
--celltype_column region_nt \
--output /data2st1/junyi/output/atac0526/dar \
--method memento-binary \
--region HIP

python 11darbysubset.py --input output/atac0416/3REGIONS_peak.h5ads \
--celltype_column region_nt \
--output /data2st1/junyi/output/atac0526/dar \
--method memento-binary \
--region AMY

