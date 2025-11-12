#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate allcools
python 04harmony.py --input output/atac1112/HIP_atacsc_neuron.h5ad --output output/atac1112 --nclust 5 --lamb 0.1
python 04harmony.py --input output/atac1112/HIP_atacsc_neuron.h5ad --output output/atac1112 --nclust 8 --lamb 0.1
python 04harmony.py --input output/atac1112/HIP_atacsc_neuron.h5ad --output output/atac1112 --nclust 10 --lamb 0.1
python 04harmony.py --input output/atac1112/HIP_atacsc_neuron.h5ad --output output/atac1112 --nclust 12 --lamb 0.1
python 04harmony.py --input output/atac1112/HIP_atacsc_neuron.h5ad --output output/atac1112 --nclust 15 --lamb 0.1

python 04harmony.py --input output/atac1112/PFC_atacsc_neuron.h5ad --output output/atac1112 --nclust 5 --lamb 0.1
python 04harmony.py --input output/atac1112/PFC_atacsc_neuron.h5ad --output output/atac1112 --nclust 8 --lamb 0.1
python 04harmony.py --input output/atac1112/PFC_atacsc_neuron.h5ad --output output/atac1112 --nclust 10 --lamb 0.1
python 04harmony.py --input output/atac1112/PFC_atacsc_neuron.h5ad --output output/atac1112 --nclust 12 --lamb 0.1
python 04harmony.py --input output/atac1112/PFC_atacsc_neuron.h5ad --output output/atac1112 --nclust 15 --lamb 0.1

python 04harmony.py --input output/atac1112/AMY_atacsc_neuron.h5ad --output output/atac1112 --nclust 5 --lamb 0.1
python 04harmony.py --input output/atac1112/AMY_atacsc_neuron.h5ad --output output/atac1112 --nclust 8 --lamb 0.1
python 04harmony.py --input output/atac1112/AMY_atacsc_neuron.h5ad --output output/atac1112 --nclust 10 --lamb 0.1
python 04harmony.py --input output/atac1112/AMY_atacsc_neuron.h5ad --output output/atac1112 --nclust 12 --lamb 0.1
python 04harmony.py --input output/atac1112/AMY_atacsc_neuron.h5ad --output output/atac1112 --nclust 15 --lamb 0.1

python 04harmony.py --input output/atac1112/HIP_atacsc_NN.h5ad --output output/atac1112 --nclust 5 --lamb 0.1
python 04harmony.py --input output/atac1112/HIP_atacsc_NN.h5ad --output output/atac1112 --nclust 10 --lamb 0.1
python 04harmony.py --input output/atac1112/HIP_atacsc_NN.h5ad --output output/atac1112 --nclust 15 --lamb 0.1

python 04harmony.py --input output/atac1112/PFC_atacsc_NN.h5ad --output output/atac1112 --nclust 5 --lamb 0.1
python 04harmony.py --input output/atac1112/PFC_atacsc_NN.h5ad --output output/atac1112 --nclust 10 --lamb 0.1
python 04harmony.py --input output/atac1112/PFC_atacsc_NN.h5ad --output output/atac1112 --nclust 15 --lamb 0.1

python 04harmony.py --input output/atac1112/AMY_atacsc_NN.h5ad --output output/atac1112 --nclust 5 --lamb 0.1
python 04harmony.py --input output/atac1112/AMY_atacsc_NN.h5ad --output output/atac1112 --nclust 10 --lamb 0.1
python 04harmony.py --input output/atac1112/AMY_atacsc_NN.h5ad --output output/atac1112 --nclust 15 --lamb 0.1