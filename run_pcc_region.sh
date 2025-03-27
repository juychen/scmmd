#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate snapatac2
cd /home/junyichen/code/scmmd
python pcc_region.py --celltype Neuron --region PFC 
python pcc_region.py --celltype Neuron --region AMY 
python pcc_region.py --celltype Sampled --region ALL 
