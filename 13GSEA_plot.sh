#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate snapatac2
# Rscript 13GSEA.r

/
# Step 2: Run Python merge and predict on generated GO results
#!/bin/bash
go_group_path="/data1st2/hannan_25/data/RNAexp/snRNA/GSEA/from_MK/2_go_term_mapping_group_20250923.xlsx"
df_class_path="/data1st2/hannan_25/data/RNAexp/snRNA/GSEA/from_MK/name_form_new.xlsx"
python_script="/home/junyichen/code/scmmd/process_GO.py"
dh="filter"
input_col='GO'
output_col='Group'
neurotransmitter_col='Neurotransmitter'
nlogp_threshold=0

# 遍历 12 个“目录”：.../*/adata_*/
for go_dir in //data2st1/junyi/output/atac1112/dar/gsea/MASTNG_dar_annotated_avg/MASTNG_dar_annotated_avg/; do
  # 输出放在该目录下的 go_gsea_heatmap/
  out_dir="${go_dir%/}/go_gsea_heatmap"
  mkdir -p "$out_dir"

  # 输出文件名用当前 adata 目录名，避免重名
  out_file="${out_dir}/$(basename "${go_dir%/}")_merged_up.csv"

  echo ">>> 处理目录: $go_dir"
  echo "    输出文件: $out_file"

  python "$python_script" \
    --go_dir_path "$go_dir" \
    --go_group_path "$go_group_path" \
    --df_class_path "$df_class_path" \
    --output_path "$out_file" \
    --duplicate_handling "$dh" \
    --input_col "$input_col" \
    --output_col "$output_col" \
    --neurotransmitter_col "$neurotransmitter_col" \
    --nlogp_threshold "$nlogp_threshold"
done
###step3 GO heatmap plot#############
/home/junyichen/anaconda3/envs/r_env/bin/Rscript /home/junyichen/code/scmmd/plotheatGSEA.r