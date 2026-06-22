#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate tobias
ulimit -n 65534

cd /data2st2/junyi/output/sn0615
# ============================================
# snRNA-seq BAM merge by condition (CON / SUS)
# Brain region: PFC (Prefrontal Cortex)
# ============================================
# snRNA BAM list: first column = condition, second column = BAM path
snrna_bams=(
"CON	/data1st2/mark/snRNA/snRNA_CUMS/cellranger_output_raw/snRNA_batch_novogene_out/cellranger_output/MW26B_PFC_novogene/outs/possorted_genome_bam.bam"
"CON	/data1st2/mark/snRNA/snRNA_CUMS/cellranger_output_raw/snRNA_batch_yunzhun_out/cellranger_output/MW26E_PFC_yunzhun/outs/possorted_genome_bam.bam"
"CON	/data1st2/mark/snRNA/snRNA_CUMS/cellranger_output_raw/snRNA_batch_beirui_out/cellranger_output/MW34C_PFC_beirui/outs/possorted_genome_bam.bam"
"CON	/data1st2/mark/snRNA/snRNA_CUMS/cellranger_output_raw/snRNA_batch_yunzhun_out/cellranger_output/MW51A_PFC_yunzhun/outs/possorted_genome_bam.bam"
"SUS	/data1st2/mark/snRNA/snRNA_CUMS/cellranger_output_raw/snRNA_batch_beirui_out/cellranger_output/MC21D_PFC_beirui/outs/possorted_genome_bam.bam"
"SUS	/data1st2/mark/snRNA/snRNA_CUMS/cellranger_output_raw/snRNA_batch_yunzhun_out/cellranger_output/MC25B_PFC_yunzhun/outs/possorted_genome_bam.bam"
"SUS	/data1st2/mark/snRNA/snRNA_CUMS/cellranger_output_raw/snRNA_batch_novogene_out/cellranger_output/MC33A_PFC_novogene/outs/possorted_genome_bam.bam"
"SUS	/data1st2/mark/snRNA/snRNA_CUMS/cellranger_output_raw/snRNA_batch_yunzhun_out/cellranger_output/MC52E_PFC_yunzhun/outs/possorted_genome_bam.bam"
)

# # Define the output folder for merged files
output_dir="/data2st2/junyi/output/sn0615/BULK_PFC"
mkdir -p "$output_dir"
mkdir -p "$output_dir/finish_list"

# # Collect BAM files by condition
# con_bams=()
# sus_bams=()
# for entry in "${snrna_bams[@]}"; do
#     cond=$(echo "$entry" | cut -f1)
#     bam=$(echo "$entry" | cut -f2)
#     if [[ "$cond" == "CON" ]]; then
#         con_bams+=("$bam")
#     elif [[ "$cond" == "SUS" ]]; then
#         sus_bams+=("$bam")
#     fi
# done

# echo "CON samples (PFC_MW): ${#con_bams[@]} BAMs"
# echo "SUS samples (PFC_MC): ${#sus_bams[@]} BAMs"

# # Merge CON -> PFC_MW.bam
# if [[ ! -f "$output_dir/finish_list/PFC_MW.txt" ]]; then
#     echo "Merging CON BAMs into PFC_MW.bam ..."
#     samtools merge -f "$output_dir/PFC_MW.bam" "${con_bams[@]}" --threads 16
#     if [[ $? -eq 0 ]]; then
#         touch "$output_dir/finish_list/PFC_MW.txt"
#         echo "  Done: PFC_MW.bam"
#     fi
# else
#     echo "PFC_MW.bam already exists, skipping."
# fi

# # Merge SUS -> PFC_MC.bam
# if [[ ! -f "$output_dir/finish_list/PFC_MC.txt" ]]; then
#     echo "Merging SUS BAMs into PFC_MC.bam ..."
#     samtools merge -f "$output_dir/PFC_MC.bam" "${sus_bams[@]}" --threads 16
#     if [[ $? -eq 0 ]]; then
#         touch "$output_dir/finish_list/PFC_MC.txt"
#         echo "  Done: PFC_MC.bam"
#     fi
# else
#     echo "PFC_MC.bam already exists, skipping."
# fi

# ============================================
# Convert merged BAM to bigWig (CPM normalized)
# ============================================
pfc_bams=("$output_dir/PFC_MW.bam" "$output_dir/PFC_MC.bam")
for bam in "${pfc_bams[@]}"; do
    base=$(basename "$bam" .bam)
    bw="${output_dir}/${base}.bw"
    flag="${output_dir}/finish_list/${base}_bw.txt"
    if [[ ! -f "$flag" ]]; then
        echo "Converting $bam -> $bw ..."
        if [[ ! -f "${bam}.bai" ]]; then
            samtools index "$bam"
        fi
        bamCoverage \
            -b "$bam" \
            -o "$bw" \
            --normalizeUsing CPM \
            --binSize 25 \
            --numberOfProcessors 64
        if [[ $? -eq 0 ]]; then
            touch "$flag"
            echo "  Done: $bw"
        fi
    else
        echo "${base}.bw already exists, skipping."
    fi
done

echo "All done. BAMs and BWs saved to $output_dir."
