#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate scanpy

#deg_path1='/data7/mark/STG/dataset/snRNA/merge_SCH/SUS_500_1000gene/merged_mast_wilcox_degs_all_log2fc0.1/merged_mast_wilcox_degs_all_log2fc0.1.csv'
#output_path1='/data7/mark/STG/dataset/snRNA/merge_SCH/SUS_500_1000gene/merged_mast_wilcox_degs_all_log2fc0.1/GO_enrichment_results2/'
#
#deg_path2='/data7/mark/STG/dataset/snRNA/merge_SCH/SUS_500_1000gene/merged_mast_wilcox_degs_all_ngenes_log2fc0.1/merged_mast_wilcox_degs_all_ngenes_log2fc0.1.csv'
#output_path2='/data7/mark/STG/dataset/snRNA/merge_SCH/SUS_500_1000gene/merged_mast_wilcox_degs_all_ngenes_log2fc0.1/GO_enrichment_results2/'

#deg_path1='/data7/mark/STG/dataset/snRNA/merge_SCH/SUS_3v3_500_1000gene/merged_mast_wilcox_degs_all_log2fc0.1/merged_mast_wilcox_degs_all_log2fc0.1.csv'
#output_path1='/data7/mark/STG/dataset/snRNA/merge_SCH/SUS_3v3_500_1000gene/merged_mast_wilcox_degs_all_log2fc0.1/GO_enrichment_results2/'
#
#deg_path2='/data7/mark/STG/dataset/snRNA/merge_SCH/SUS_3v3_500_1000gene/merged_mast_wilcox_degs_all_ngenes_log2fc0.1/merged_mast_wilcox_degs_all_ngenes_log2fc0.1.csv'
#output_path2='/data7/mark/STG/dataset/snRNA/merge_SCH/SUS_3v3_500_1000gene/merged_mast_wilcox_degs_all_ngenes_log2fc0.1/GO_enrichment_results2/'

#deg_path1='/data7/mark/STG/dataset/snRNA/merge_SCH/RES_3v3_500_1000gene/merged_mast_wilcox_degs_log2fc0.1/merged_mast_wilcox_degs_log2fc0.1.csv'
#output_path1=$(dirname "$deg_path1")
#
#deg_path2='/data7/mark/STG/dataset/snRNA/merge_SCH/RES_3v3_500_1000gene/merged_mast_wilcox_degs_all_ngenes_log2fc0.1/merged_mast_wilcox_degs_all_ngenes_log2fc0.1.csv'
#output_path2=$(dirname "$deg_path2")
#
#regions=("MB" "HY" "TH" "AMY" "HPF" "STR" "PFC")
#for region in "${regions[@]}"; do
#  echo "$deg_path1"
#  Rscript gost_tf.r -i "$deg_path1" -o "$output_path1" -r "$region" &
#done
#wait
#
#for region in "${regions[@]}"; do
#  Rscript gost_tf.r -i "$deg_path2" -o "$output_path2" -r "$region" &
#done
#wait



## Define DEG paths
#deg_paths=(
#  "/data7/mark/STG/dataset/snRNA/merge_SCH/RES_3v3_500_1000gene/merged_mast_wilcox_degs_log2fc0.1/merged_mast_wilcox_degs_log2fc0.1.csv"
#  "/data7/mark/STG/dataset/snRNA/merge_SCH/RES_3v3_500_1000gene/merged_mast_wilcox_degs_all_ngenes_log2fc0.1/merged_mast_wilcox_degs_all_ngenes_log2fc0.1.csv"
#)
#
## Define regions
#regions=("MB" "HY" "TH" "AMY" "HPF" "STR" "PFC" "iCTX")
#go_rscript_path="/data7/mark/STG/dataset/snRNA/merge_SCH/code_go_mast/gost_tf.r"
#
## Loop over DEG paths
#for deg_path in "${deg_paths[@]}"; do
#  base_output=$(dirname "$deg_path")
#  output_path="${base_output}/GO_enrichment_results2"
#
#  # Loop over regions
#  for region in "${regions[@]}"; do
#    echo "Running gost_tf.r for $region with input $deg_path"
#    Rscript "$go_rscript_path" -i "$deg_path" -o "$output_path" -r "$region" &
#  done
#  wait   # wait for all regions of this deg_path
#done




# Define DEG paths
#deg_paths=(
#  "/data7/mark/STG/dataset/snRNA/merge_SCH/SUS_500_1000gene/merged_mast_wilcox_degs_all_fdr_log2fc0.1/merged_mast_wilcox_degs_all_fdr_log2fc0.1.csv"
#  "/data7/mark/STG/dataset/snRNA/merge_SCH/SUS_500_1000gene/merged_mast_wilcox_degs_all_ngenes_fdr_log2fc0.1/merged_mast_wilcox_degs_all_ngenes_fdr_log2fc0.1.csv"
#  "/data7/mark/STG/dataset/snRNA/merge_SCH/SUS_500_1000gene/merged_mast_wilcox_degs_all_ngenes_saturation_fdr_log2fc0.1/merged_mast_wilcox_degs_all_ngenes_saturation_fdr_log2fc0.1.csv"
#  "/data7/mark/STG/dataset/snRNA/merge_SCH/SUS_3v3_500_1000gene/merged_mast_wilcox_degs_all_fdr_log2fc0.1/merged_mast_wilcox_degs_all_fdr_log2fc0.1.csv"
#  "/data7/mark/STG/dataset/snRNA/merge_SCH/SUS_3v3_500_1000gene/merged_mast_wilcox_degs_all_ngenes_fdr_log2fc0.1/merged_mast_wilcox_degs_all_ngenes_fdr_log2fc0.1.csv"
#  "/data7/mark/STG/dataset/snRNA/merge_SCH/SUS_3v3_500_1000gene/merged_mast_wilcox_degs_all_ngenes_saturation_fdr_log2fc0.1/merged_mast_wilcox_degs_all_ngenes_saturation_fdr_log2fc0.1.csv"
#  "/data7/mark/STG/dataset/snRNA/merge_SCH/RES_3v3_500_1000gene/merged_mast_wilcox_degs_all_fdr_log2fc0.1/merged_mast_wilcox_degs_all_fdr_log2fc0.1.csv"
#  "/data7/mark/STG/dataset/snRNA/merge_SCH/RES_3v3_500_1000gene/merged_mast_wilcox_degs_all_ngenes_fdr_log2fc0.1/merged_mast_wilcox_degs_all_ngenes_fdr_log2fc0.1.csv"
#  "/data7/mark/STG/dataset/snRNA/merge_SCH/RES_3v3_500_1000gene/merged_mast_wilcox_degs_all_ngenes_saturation_fdr_log2fc0.1/merged_mast_wilcox_degs_all_ngenes_saturation_fdr_log2fc0.1.csv"
#  "/data7/mark/STG/dataset/snRNA/merge_SCH/RES_3v3_beirui_500_1000gene/merged_mast_wilcox_degs_beirui_fdr_log2fc0.1/merged_mast_wilcox_degs_beirui_fdr_log2fc0.1.csv"
#  "/data7/mark/STG/dataset/snRNA/merge_SCH/RES_3v3_beirui_500_1000gene/merged_mast_wilcox_degs_beirui_ngenes_fdr_log2fc0.1/merged_mast_wilcox_degs_beirui_ngenes_fdr_log2fc0.1.csv"
#  "/data7/mark/STG/dataset/snRNA/merge_SCH/RES_3v3_beirui_500_1000gene/merged_mast_wilcox_degs_beirui_ngenes_saturation_fdr_log2fc0.1/merged_mast_wilcox_degs_beirui_ngenes_saturation_fdr_log2fc0.1.csv"
#  "/data7/mark/STG/dataset/snRNA/merge_SCH/CSDS_3v3_500_1000gene/merged_mast_wilcox_degs_beirui_fdr_log2fc0.1/merged_mast_wilcox_degs_beirui_fdr_log2fc0.1.csv"
#  "/data7/mark/STG/dataset/snRNA/merge_SCH/CSDS_3v3_500_1000gene/merged_mast_wilcox_degs_beirui_ngenes_fdr_log2fc0.1/merged_mast_wilcox_degs_beirui_ngenes_fdr_log2fc0.1.csv"
#  "/data7/mark/STG/dataset/snRNA/merge_SCH/CSDS_3v3_500_1000gene/merged_mast_wilcox_degs_beirui_ngenes_saturation_fdr_log2fc0.1/merged_mast_wilcox_degs_beirui_ngenes_saturation_fdr_log2fc0.1.csv"
#  "/data7/mark/STG/dataset/snRNA/merge_SCH/CSDS_3v3_500_1000gene/merged_mast_wilcox_degs_All_log2fc0.1/merged_mast_wilcox_degs_All_log2fc0.1.csv"
#  "/data7/mark/STG/dataset/snRNA/merge_SCH/CSDS_3v3_500_1000gene/merged_mast_wilcox_degs_Allng_log2fc0.1/merged_mast_wilcox_degs_Allng_log2fc0.1.csv"
#  "/data7/mark/STG/dataset/snRNA/merge_SCH/CSDS_3v3_500_1000gene/merged_mast_wilcox_degs_Allsa_log2fc0.1/merged_mast_wilcox_degs_Allsa_log2fc0.1.csv"
#"/data7/mark/STG/dataset/snRNA/merge_SCH/CSDS_3v3_500_1000gene/merged_mast_wilcox_degs_all_fdr_log2fc0.1/merged_mast_wilcox_degs_all_fdr_log2fc0.1.csv"
#"/data7/mark/STG/dataset/snRNA/merge_SCH/CSDS_3v3_500_1000gene/merged_mast_wilcox_degs_all_ngenes_fdr_log2fc0.1/merged_mast_wilcox_degs_all_ngenes_fdr_log2fc0.1.csv"
#"/data7/mark/STG/dataset/snRNA/merge_SCH/CSDS_4v3_500_1000gene/merged_mast_wilcox_degs_all_fdr_log2fc0.1/merged_mast_wilcox_degs_all_fdr_log2fc0.1.csv"
#"/data7/mark/STG/dataset/snRNA/merge_SCH/CSDS_4v3_500_1000gene/merged_mast_wilcox_degs_all_ngenes_fdr_log2fc0.1/merged_mast_wilcox_degs_all_ngenes_fdr_log2fc0.1.csv"
#"/data7/mark/STG/dataset/snRNA/merge_SCH/RES_3v3_500_1000gene_new/merged_mast_wilcox_degs_all_fdr_log2fc0.1/merged_mast_wilcox_degs_all_fdr_log2fc0.1.csv"
#"/data7/mark/STG/dataset/snRNA/merge_SCH/RES_3v3_500_1000gene_new/merged_mast_wilcox_degs_all_ngenes_fdr_log2fc0.1/merged_mast_wilcox_degs_all_ngenes_fdr_log2fc0.1.csv"
#"/data7/mark/STG/dataset/snRNA/merge_SCH/RES_4v3_500_1000gene/merged_mast_wilcox_degs_all_fdr_log2fc0.1/merged_mast_wilcox_degs_all_fdr_log2fc0.1.csv"
#"/data7/mark/STG/dataset/snRNA/merge_SCH/RES_4v3_500_1000gene/merged_mast_wilcox_degs_all_ngenes_fdr_log2fc0.1/merged_mast_wilcox_degs_all_ngenes_fdr_log2fc0.1.csv"
#'/data7/mark/STG/dataset/snRNA/merge_SCH/SUS_500_1000gene/merged_mast_wilcox_degs_all_ngenes_only_fdr_log2fc0.1/merged_mast_wilcox_degs_all_ngenes_only_fdr_log2fc0.1.csv'
#'/data7/mark/STG/dataset/snRNA/merge_SCH/SUS_3v3_500_1000gene/merged_mast_wilcox_degs_all_ngenes_only_fdr_log2fc0.1/merged_mast_wilcox_degs_all_ngenes_only_fdr_log2fc0.1.csv'
#'/data7/mark/STG/dataset/snRNA/merge_SCH/CSDS_4v3_500_1000gene/merged_memento_wilcox_mast_degs_all_tech_log2fc0.2/merged_memento_wilcox_mast_degs_all_tech_log2fc0.2.csv'
#'/data7/mark/STG/dataset/snRNA/merge_SCH/CSDS_3v3_500_1000gene/merged_memento_wilcox_mast_degs_all_tech_log2fc0.2/merged_memento_wilcox_mast_degs_all_tech_log2fc0.2.csv'
#'/data7/mark/STG/dataset/snRNA/merge_SCH/CSDS_3v3_beirui_500_1000gene/merged_memento_wilcox_mast_degs_all_tech_log2fc0.2/merged_memento_wilcox_mast_degs_all_tech_log2fc0.2.csv'
#)

# base_dirs=(
#   "/data7/mark/STG/dataset/snRNA/merge_SCH_new/SUS_4v4_500_1000gene"
#   # "/data7/mark/STG/dataset/snRNA/merge_SCH_new/SUS_3v3_500_1000gene"
#   "/data7/mark/STG/dataset/snRNA/merge_SCH_new/RES_4v3_500_1000gene"
#   # "/data7/mark/STG/dataset/snRNA/merge_SCH_new/RES_3v3_500_1000gene"
#   "/data7/mark/STG/dataset/snRNA/merge_SCH_new/RES_3v3_500_1000gene_beirui"
#   "/data7/mark/STG/dataset/snRNA/merge_SCH_new/CSDS_4v3_500_1000gene"
#   # "/data7/mark/STG/dataset/snRNA/merge_SCH_new/CSDS_3v3_500_1000gene"
#   "/data7/mark/STG/dataset/snRNA/merge_SCH_new/CSDS_3v3_500_1000gene_beirui"
#  "/data7/mark/STG/dataset/snRNA/merge_SCH_new/CSRES_4v3_500_1000gene"
# #  "/data7/mark/STG/dataset/snRNA/merge_SCH_new/CSRES_3v3_500_1000gene"
#  "/data7/mark/STG/dataset/snRNA/merge_SCH_new/CSRES_3v3_500_1000gene_beirui"
# )

# base_dirs=(
#   # "/data11/junyi/stg1028/CSDS_3VB/"
#   # "/data11/junyi/stg1028/CSDS_4VN/"
#   # "/data11/junyi/stg1028/CSRES_3VB/"
#   # "/data11/junyi/stg1028/CSRES_4VN/"
#   "/data2st2/junyi/output/stg1028/CUMS_4VN/"
#   # "/data11/junyi/stg1028/CURES_3VB/"
#   # "/data11/junyi/stg1028/CURES_4VN/"
# )

# #deg_dirs=(
# #	'merged_memento_mast3_cov4_ratio_degs_bonferroni_all_company_ngeneson_avg_pct_counts_mt_avg_pct_counts_ribo_avg_sample_downsampled_ratio_log2fc0.1'
# #	'merged_memento_mast3_cov4_ratio_degs_fdr_all_company_ngeneson_avg_pct_counts_mt_avg_pct_counts_ribo_avg_sample_downsampled_ratio_log2fc0.1'
# #	'merged_memento_mast3_cov4_ratio_wilcox_degs_bonferroni_all_company_ngeneson_avg_pct_counts_mt_avg_pct_counts_ribo_avg_sample_downsampled_ratio_log2fc0.1'
# #	'merged_memento_mast3_cov4_ratio_wilcox_degs_fdr_all_company_ngeneson_avg_pct_counts_mt_avg_pct_counts_ribo_avg_sample_downsampled_ratio_log2fc0.1'          
# #)
# deg_dirs=(
# #	'merged_mast_wilcox_ngsamtrb_degs_fdr_log2fc0.1'
# #	'merged_mast_wilcox_ngsamtrb_degs_bonferroni_log2fc0.1'
# 	# 'merged_mast_ext_wilcox_ngsamtrb_ext_degs_fdr_log2fc0.1'
# 	# 'merged_mast_ext_wilcox_ngsamtrb_ext_degs_bonferroni_log2fc0.1'
# #	'merged_memento_mast3_cov4_ratio_wilcox_degs_bonferroni_all_company_ngeneson_avg_pct_counts_mt_avg_pct_counts_ribo_avg_sample_downsampled_ratio_log2fc0.1'
# #	'merged_memento_mast3_cov4_ratio_wilcox_degs_fdr_all_company_ngeneson_avg_pct_counts_mt_avg_pct_counts_ribo_avg_sample_downsampled_ratio_log2fc0.1'          
# "mast_ngsa_degs_fdr_log2fc0.1"
# "mast_ngsa_degs_fdr_log2fc0.15"
# "mast_ngsa_degs_fdr_log2fc0.2"
# "mast_ngsarb_degs_fdr_log2fc0.1"
# "mast_ngsarb_degs_fdr_log2fc0.15"
# "mast_ngsarb_degs_fdr_log2fc0.2"
# )

# csv_suffix="filtered"
# deg_paths=()
# for base in "${base_dirs[@]}"; do
#   for deg in "${deg_dirs[@]}"; do
    
#     # If the base path does NOT contain "SUS", remove "company" from deg
#     if [[ "$base" != *"CUMS"* ]]; then
#       deg_mod="${deg//_company/}"
#     else
#       deg_mod="$deg"
#     fi

#     # Primary candidate file
#     candidate="$base/$deg_mod/$deg_mod.csv"

#     if [[ -f "$candidate" ]]; then
#       deg_paths+=("$candidate")
#     else
#       echo "Warning: DEG file not found -> $candidate"
#     fi

#     # If csv_suffix is not empty, add an additional candidate file
#     if [[ -n "$csv_suffix" ]]; then
#       candidate_suffix="$base/$deg_mod/${deg_mod}_${csv_suffix}.csv"
#       if [[ -f "$candidate_suffix" ]]; then
#         deg_paths+=("$candidate_suffix")
#       else
#         echo "Warning: DEG file not found -> $candidate_suffix"
#       fi
#     fi

#   done
# done


# Define regions
regions=("AMY" "HPF" "PFC")

# Path to R script
go_rscript_path="/home/junyichen/code/scmmd/custerprofiler.r"

# Python script
python_script="/data2st2/junyi/code/sn/example_go_merge_and_predict3.py"

# Mapping and class info
go_group_path="/data2st2/junyi/code/sn/2_go_term_mapping_group_20250923.xlsx"
df_class_path="/data2st2/junyi/code/sn/name_form_new.xlsx"

# Duplicate handling strategies
duplicate_handlings=("filter" "Down" "Up" "nlogp_keep_highest" "nlogp_compare")
input_col='GO'
output_col='Group'
neurotransmitter_col='Neurotransmitter'
nlogp_threshold=0.1

# Max parallel jobs
max_jobs=8
job_count=0

deg_paths=("/data2st1/junyi/output/atac1112/dar/celltype.L2/MASTNG_dar_annotated_avg.csv")
# Step 1: Run R script to generate GO enrichment results
for deg_path in "${deg_paths[@]}"; do
    parent_dir=$(dirname "$deg_path")
    deg_filename=$(basename "$deg_path" .csv)
    base_output="${parent_dir}/${deg_filename}"
    mkdir -p "$base_output"
    
    go_output_dir="${base_output}/CP_enrichment_results2"

#    if [[ -d "$go_output_dir" && -n "$(ls -A "$go_output_dir"/*.csv 2>/dev/null)" ]]; then
#        echo "Skipping $go_output_dir (CSV files already exist)"
#        continue
#    fi

    mkdir -p "$go_output_dir"

    conda activate r_env
    echo "Activated r_env for R script execution, output dir: $go_output_dir"

    for region in "${regions[@]}"; do
    
        echo "Running gost_tf.r for $region with input $deg_path"
        Rscript "$go_rscript_path" -i "$deg_path" -o "$go_output_dir" -r "$region" &

        ((job_count++))
        if ((job_count % max_jobs == 0)); then
            wait
        fi
    done
done
wait
echo "=== All GO enrichment results generated ==="

# Step 2: Run Python merge and predict on generated GO results
for deg_path in "${deg_paths[@]}"; do
    parent_dir=$(dirname "$deg_path")
    deg_filename=$(basename "$deg_path" .csv)
    base_output="${parent_dir}/${deg_filename}"
    go_dir="${base_output}/CP_enrichment_results2"

    for dh in "${duplicate_handlings[@]}"; do
        output_dir="${base_output}/${dh}"
        mkdir -p "$output_dir"
        output_file="${output_dir}/2_go_merged_with_predictions.csv"

        echo "Running Python merge/predict for GO results in $go_dir with duplicate_handling=$dh"
        python "$python_script" \
            --go_dir_path "$go_dir" \
            --go_group_path "$go_group_path" \
            --df_class_path "$df_class_path" \
            --output_path "$output_file" \
            --duplicate_handling "$dh" \
            --input_col "$input_col" \
            --output_col "$output_col" \
            --neurotransmitter_col "$neurotransmitter_col" \
            --nlogp_threshold "$nlogp_threshold"
    done
done

echo "=== All merge and prediction tasks completed ==="



