#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh

conda activate snapatac2

# output path make if not exists
mkdir -p /data2st1/junyi/output/atac0627/darmr

python 11darbysubset.py --input /data2st1/junyi/output/atac0627/3REGIONS_dmr.h5ads \
--output /data2st1/junyi/output/atac0627/darmr \
--celltype_column All \
--method memento-binary \
--region PFC

python 11darbysubset.py --input /data2st1/junyi/output/atac0627/3REGIONS_dmr.h5ads \
--celltype_column All \
--output /data2st1/junyi/output/atac0627/darmr \
--method memento-binary \
--region HIP

python 11darbysubset.py --input /data2st1/junyi/output/atac0627/3REGIONS_dmr.h5ads \
--celltype_column All \
--output /data2st1/junyi/output/atac0627/dar \
--method memento-binary \
--region AMY

python 11darbysubset.py --input /data2st1/junyi/output/atac0627/3REGIONS_peak.h5ads \
--output /data2st1/junyi/output/atac0627/dar \
--celltype_column celltype.L2 \
--method wilcoxon \
--region PFC

python 11darbysubset.py --input /data2st1/junyi/output/atac0627/3REGIONS_peak.h5ads \
--celltype_column celltype.L2 \
--output /data2st1/junyi/output/atac0627/dar \
--method wilcoxon \
--region HIP

python 11darbysubset.py --input /data2st1/junyi/output/atac0627/3REGIONS_peak.h5ads \
--celltype_column celltype.L2 \
--output /data2st1/junyi/output/atac0627/dar \
--method wilcoxon \
--region AMY
