#!/bin/bash
# python harmony.py --nclust 8 --lamb 0.1
# python harmony.py --nclust 10 --lamb 0.1
# python harmony.py --nclust 12 --lamb 0.1
# python harmony.py --nclust 8 --lamb 0.5
# python harmony.py --nclust 10 --lamb 0.5
# python harmony.py --nclust 12 --lamb 0.5
# python harmony.py --nclust 8 --lamb 0.3
# python harmony.py --nclust 10 --lamb 0.3
# python harmony.py --nclust 12 --lamb 0.3
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate snapatac2
cd /home/junyichen/code/scmmd
python 07iterativeclust.py
