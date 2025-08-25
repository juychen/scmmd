#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate allcools

python 04harmony.py --input /data1st2/junyi/output/methlypub0823/PFC_atacsc_neuron.h5ad --output /data1st2/junyi/output/methlypub0823 --nclust 15 --lamb 0.1
# python 04harmony.py --input /data1st2/junyi/output/methlypub0823/AMY_atacsc_neuron.h5ad --output /data1st2/junyi/output/methlypub0823 --nclust 15 --lamb 0.1
# python 04harmony.py --input /data1st2/junyi/output/methlypub0823/HPF_atacsc_neuron.h5ad --output /data1st2/junyi/output/methlypub0823 --nclust 15 --lamb 0.1


# python 04harmony.py --input /data1st2/junyi/output/methlypub0823/PFC_atacsc_NN.h5ad --output /data1st2/junyi/output/methlypub0823 --nclust 15 --lamb 0.1
# python 04harmony.py --input /data1st2/junyi/output/methlypub0823/AMY_atacsc_NN.h5ad --output /data1st2/junyi/output/methlypub0823 --nclust 15 --lamb 0.1
# python 04harmony.py --input /data1st2/junyi/output/methlypub0823/HPF_atacsc_NN.h5ad --output /data1st2/junyi/output/methlypub0823 --nclust 15 --lamb 0.1


