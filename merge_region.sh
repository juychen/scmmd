#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh

conda activate snapatac2
bedtools merge -i /data2st1/junyi/snmc_bedgraph/merged/all_unified.bedGraph -d 150 -c 4,5,6,...,277 -o mean > /data2st1/junyi/snmc_bedgraph/merged/merged_output.bedGraph
