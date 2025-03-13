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
python harmony.py --input atacsc-3region-clustering.h5ad --nclust 5 --lamb 0.1
python harmony.py --input atacsc-3region-clustering.h5ad --nclust 8 --lamb 0.1
python harmony.py --input atacsc-3region-clustering.h5ad --nclust 10 --lamb 0.1
python harmony.py --input atacsc-3region-clustering.h5ad --nclust 12 --lamb 0.1
python harmony.py --input atacsc-3region-clustering.h5ad --nclust 15 --lamb 0.31
