#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate allcools
# python 04harmony.py --input output/atac0627/HIP_atlassc_NN.h5ad --output output/atac0627 --nclust 5 --lamb 0.1
# python 04harmony.py --input output/atac0627/HIP_atlassc_NN.h5ad --output output/atac0627 --nclust 8 --lamb 0.1
python 04harmony.py --input output/atac0627/HIP_atlassc_NN.h5ad --output output/atac0627 --nclust 10 --lamb 0.1
# python 04harmony.py --input output/atac0627/HIP_atlassc_NN.h5ad --output output/atac0627 --nclust 12 --lamb 0.1
python 04harmony.py --input output/atac0627/HIP_atlassc_NN.h5ad --output output/atac0627 --nclust 15 --lamb 0.1

# python 04harmony.py --input output/atac0627/PFC_atlassc_NN.h5ad --output output/atac0627 --nclust 5 --lamb 0.1
# python 04harmony.py --input output/atac0627/PFC_atlassc_NN.h5ad --output output/atac0627 --nclust 8 --lamb 0.1
python 04harmony.py --input output/atac0627/PFC_atlassc_NN.h5ad --output output/atac0627 --nclust 10 --lamb 0.1
# python 04harmony.py --input output/atac0627/PFC_atlassc_NN.h5ad --output output/atac0627 --nclust 12 --lamb 0.1
python 04harmony.py --input output/atac0627/PFC_atlassc_NN.h5ad --output output/atac0627 --nclust 15 --lamb 0.1

# python 04harmony.py --input output/atac0627/AMY_atlassc_NN.h5ad --output output/atac0627 --nclust 5 --lamb 0.1
# python 04harmony.py --input output/atac0627/AMY_atlassc_NN.h5ad --output output/atac0627 --nclust 8 --lamb 0.1
python 04harmony.py --input output/atac0627/AMY_atlassc_NN.h5ad --output output/atac0627 --nclust 10 --lamb 0.1
# python 04harmony.py --input output/atac0627/AMY_atlassc_NN.h5ad --output output/atac0627 --nclust 12 --lamb 0.1
python 04harmony.py --input output/atac0627/AMY_atlassc_NN.h5ad --output output/atac0627 --nclust 15 --lamb 0.1

# python 04harmony.py --input output/atac0627/HIP_atlassc_neuron.h5ad --output output/atac0627 --nclust 5 --lamb 0.1
# python 04harmony.py --input output/atac0627/HIP_atlassc_neuron.h5ad --output output/atac0627 --nclust 8 --lamb 0.1
python 04harmony.py --input output/atac0627/HIP_atlassc_neuron.h5ad --output output/atac0627 --nclust 10 --lamb 0.1
# python 04harmony.py --input output/atac0627/HIP_atlassc_neuron.h5ad --output output/atac0627 --nclust 12 --lamb 0.1
python 04harmony.py --input output/atac0627/HIP_atlassc_neuron.h5ad --output output/atac0627 --nclust 15 --lamb 0.1

# python 04harmony.py --input output/atac0627/PFC_atlassc_neuron.h5ad --output output/atac0627 --nclust 5 --lamb 0.1
# python 04harmony.py --input output/atac0627/PFC_atlassc_neuron.h5ad --output output/atac0627 --nclust 8 --lamb 0.1
python 04harmony.py --input output/atac0627/PFC_atlassc_neuron.h5ad --output output/atac0627 --nclust 10 --lamb 0.1
# python 04harmony.py --input output/atac0627/PFC_atlassc_neuron.h5ad --output output/atac0627 --nclust 12 --lamb 0.1
python 04harmony.py --input output/atac0627/PFC_atlassc_neuron.h5ad --output output/atac0627 --nclust 15 --lamb 0.1

# python 04harmony.py --input output/atac0627/AMY_atlassc_neuron.h5ad --output output/atac0627 --nclust 5 --lamb 0.1
# python 04harmony.py --input output/atac0627/AMY_atlassc_neuron.h5ad --output output/atac0627 --nclust 8 --lamb 0.1
python 04harmony.py --input output/atac0627/AMY_atlassc_neuron.h5ad --output output/atac0627 --nclust 10 --lamb 0.1
# python 04harmony.py --input output/atac0627/AMY_atlassc_neuron.h5ad --output output/atac0627 --nclust 12 --lamb 0.1
python 04harmony.py --input output/atac0627/AMY_atlassc_neuron.h5ad --output output/atac0627 --nclust 15 --lamb 0.1


