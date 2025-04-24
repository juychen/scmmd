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
python harmony.py --input output/atac0416/AMY_atacsc_neuron.h5ad --output output/atac0416 --nclust 5 --lamb 0.1
python harmony.py --input output/atac0416/AMY_atacsc_neuron.h5ad --output output/atac0416 --nclust 8 --lamb 0.1
python harmony.py --input output/atac0416/AMY_atacsc_neuron.h5ad --output output/atac0416 --nclust 10 --lamb 0.1
python harmony.py --input output/atac0416/AMY_atacsc_neuron.h5ad --output output/atac0416 --nclust 12 --lamb 0.1
python harmony.py --input output/atac0416/AMY_atacsc_neuron.h5ad --output output/atac0416 --nclust 15 --lamb 0.31

python harmony.py --input output/atac0416/PFC_atacsc_neuron.h5ad --output output/atac0416 --nclust 5 --lamb 0.1
python harmony.py --input output/atac0416/PFC_atacsc_neuron.h5ad --output output/atac0416 --nclust 8 --lamb 0.1
python harmony.py --input output/atac0416/PFC_atacsc_neuron.h5ad --output output/atac0416 --nclust 10 --lamb 0.1
python harmony.py --input output/atac0416/PFC_atacsc_neuron.h5ad --output output/atac0416 --nclust 12 --lamb 0.1
python harmony.py --input output/atac0416/PFC_atacsc_neuron.h5ad --output output/atac0416 --nclust 15 --lamb 0.31
