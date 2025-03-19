#!/bin/bash
ts=$(date +"%Y-%m-%d-%H-%M-%S")
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate snapatac2
python snap_dar.py