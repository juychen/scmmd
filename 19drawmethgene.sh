#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate r_env
# Rscript gost_tfmodule.r 'AMY'
# regions = ['HPF','Isocortex', 'AMY',
#        'PFC', 'TH', 'STR',  'HY', 'MB']
Rscript gost_tfmodule.r 'NN' &\
Rscript gost_tfmodule.r 'TH' &\
Rscript gost_tfmodule.r 'HPF' &\
Rscript gost_tfmodule.r 'Isocortex' ;\
Rscript gost_tfmodule.r 'STR' &\
Rscript gost_tfmodule.r 'HY' &\
Rscript gost_tfmodule.r 'MB' &\
Rscript gost_tfmodule.r 'PFC'
