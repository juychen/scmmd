#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate r_env
#Rscript gost_tf.r 'AMY'
# regions = ['HPF','Isocortex', 'AMY',
#        'PFC', 'TH', 'STR',  'HY', 'MB']
Rscript gost_tf.r 'NN'

# Rscript gost_tf.r 'TH' &\
# Rscript gost_tf.r 'HPF' &\
# Rscript gost_tf.r 'Isocortex' &\
# Rscript gost_tf.r 'STR' ;\
# Rscript gost_tf.r 'HY' &\
# Rscript gost_tf.r 'MB' &\
# Rscript gost_tf.r 'PFC'
