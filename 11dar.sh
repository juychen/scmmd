#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh

conda activate snapatac2

python 11darbysubset.py --input output/atac0627/3REGIONS_peak.h5ads \
--output /data2st1/junyi/output/atac0627/dar \
--celltype_column celltype.L2 \
--method memento-binary \
--region PFC

python 11darbysubset.py --input output/atac0627/3REGIONS_peak.h5ads \
--celltype_column celltype.L2 \
--output /data2st1/junyi/output/atac0627/dar \
--method memento-binary \
--region HIP

python 11darbysubset.py --input output/atac0627/3REGIONS_peak.h5ads \
--celltype_column celltype.L2 \
--output /data2st1/junyi/output/atac0627/dar \
--method memento-binary \
--region AMY

