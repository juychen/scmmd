#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate scvi-env
output=/data2st2/junyi/output/stg1028/scvisn/L2

# python 11subsetloom_new.py --input TH=/data3/junyi/scvi/TH.h5ad --output "$output" --layer "scvi_reconstructed_counts" &
# python 11subsetloom_new.py --input STR=/data3/junyi/scvi/STR.h5ad --output "$output" --layer "scvi_reconstructed_counts" &
# python 11subsetloom_new.py --input PFC=/data3/junyi/scvi/PFC.h5ad --output "$output" --layer "scvi_reconstructed_counts" &
# python 11subsetloom_new.py --input MB=/data3/junyi/scvi/MB.h5ad --output "$output" --layer "scvi_reconstructed_counts" &
# wait

# python 11subsetloom_new.py --input iCTX=/data3/junyi/scvi/iCTX.h5ad --output "$output" --layer "scvi_reconstructed_counts" &
# python 11subsetloom_new.py --input HY=/data3/junyi/scvi/HY.h5ad --output "$output" --layer "scvi_reconstructed_counts" &
# python 11subsetloom_new.py --input HPF=/data3/junyi/scvi/HPF.h5ad --output "$output" --layer "scvi_reconstructed_counts" &
python 11subsetloom_new.py --input AMY=/data3/junyi/scvi/AMY.h5ad --output "$output" --layer "scvi_reconstructed_counts" &
wait