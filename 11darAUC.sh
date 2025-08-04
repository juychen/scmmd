#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh

conda activate snapatac2

# output path make if not exists
mkdir -p /data2st1/junyi/output/atac0627/darauc

python 11darbysubset.py --input /data2st1/junyi/output/atac0627/3REGIONS_dmr.h5ads \
--output /data2st1/junyi/output/atac0627/darauc \
--celltype_column All \
--method memento-binary \
--region PFC