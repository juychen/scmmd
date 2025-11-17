#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate scenicplus
ulimit -n 65534
python ctx.py