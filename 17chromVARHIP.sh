#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate scenicplus
ulimit -n 65534

script_path="/home/junyichen/code/scmmd/test_chromvar.py"
output_root="/data2st1/junyi/output/atac1112/chromvar/HIP"

input_files=(
    "/data1st2/junyi/output/atac1112/subset/region_nt/HIP_HIP_GABA.h5ad"
    "/data1st2/junyi/output/atac1112/subset/region_nt/HIP_HIP_Glut.h5ad"
    "/data1st2/junyi/output/atac1112/subset/region_nt/HIP_NN.h5ad"
)

mkdir -p "$output_root"

for input_h5ad in "${input_files[@]}"; do
    sample_name=$(basename "$input_h5ad" .h5ad)
    sample_outdir="$output_root/$sample_name"
    output_h5ad="$sample_outdir/${sample_name}_chromvar.h5ad"

    mkdir -p "$sample_outdir"
    echo "Running chromVAR for $sample_name"
    python "$script_path" "$input_h5ad" "$output_h5ad"
done

echo "chromVAR complete. Output saved under $output_root."
