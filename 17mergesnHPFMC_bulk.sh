#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate tobias
ulimit -n 65534

cd /data2st2/junyi/output/sn0615
# ============================================
# snRNA-seq BAM merge by condition (CON / SUS)
# Brain region: HPF (Hippocampal Formation)
# ============================================
# snRNA BAM list: first column = condition, second column = BAM path
snrna_bams=(
"CON	/data1st2/mark/snRNA/snRNA_CUMS/cellranger_output_raw/snRNA_batch_yunzhun_out/cellranger_output/MW22B_HIP_yunzhun/outs/possorted_genome_bam.bam"
"CON	/data1st2/mark/snRNA/snRNA_CUMS/cellranger_output_raw/snRNA_batch_yunzhun_out/cellranger_output/MW45A_HIP_yunzhun/outs/possorted_genome_bam.bam"
"CON	/data1st2/mark/snRNA/snRNA_CUMS/cellranger_output_raw/snRNA_batch_beirui_out/cellranger_output/MW45C_HIP_beirui/outs/possorted_genome_bam.bam"
"CON	/data1st2/mark/snRNA/snRNA_CUMS/cellranger_output_raw/snRNA_batch_yunzhun_out/cellranger_output/MW51A_HIP_yunzhun/outs/possorted_genome_bam.bam"
"SUS	/data1st2/mark/snRNA/snRNA_CUMS/cellranger_output_raw/snRNA_batch_yunzhun_out/cellranger_output/MC25B_HIP_yunzhun/outs/possorted_genome_bam.bam"
"SUS	/data1st2/mark/snRNA/snRNA_CUMS/cellranger_output_raw/snRNA_batch_yunzhun_out/cellranger_output/MC33A_HIP_yunzhun/outs/possorted_genome_bam.bam"
"SUS	/data1st2/mark/snRNA/snRNA_CUMS/cellranger_output_raw/snRNA_batch_yunzhun_out/cellranger_output/MC48D_HIP_yunzhun/outs/possorted_genome_bam.bam"
"SUS	/data1st2/mark/snRNA/snRNA_CUMS/cellranger_output_raw/snRNA_batch_beirui_out/cellranger_output/MC48E_HIP_beirui/outs/possorted_genome_bam.bam"
)

# Define the output folder for merged files
output_dir="/data2st2/junyi/output/sn0615/BULK_HPF"
mkdir -p "$output_dir"
mkdir -p "$output_dir/finish_list"

# Collect BAM files by condition
con_bams=()
sus_bams=()
for entry in "${snrna_bams[@]}"; do
    cond=$(echo "$entry" | cut -f1)
    bam=$(echo "$entry" | cut -f2)
    if [[ "$cond" == "CON" ]]; then
        con_bams+=("$bam")
    elif [[ "$cond" == "SUS" ]]; then
        sus_bams+=("$bam")
    fi
done

echo "CON samples (HPF_MW): ${#con_bams[@]} BAMs"
echo "SUS samples (HPF_MC): ${#sus_bams[@]} BAMs"

# Merge CON -> HPF_MW.bam
if [[ ! -f "$output_dir/finish_list/HPF_MW.txt" ]]; then
    echo "Merging CON BAMs into HPF_MW.bam ..."
    samtools merge -f "$output_dir/HPF_MW.bam" "${con_bams[@]}" --threads 16
    if [[ $? -eq 0 ]]; then
        touch "$output_dir/finish_list/HPF_MW.txt"
        echo "  Done: HPF_MW.bam"
    fi
else
    echo "HPF_MW.bam already exists, skipping."
fi

# Merge SUS -> HPF_MC.bam
if [[ ! -f "$output_dir/finish_list/HPF_MC.txt" ]]; then
    echo "Merging SUS BAMs into HPF_MC.bam ..."
    samtools merge -f "$output_dir/HPF_MC.bam" "${sus_bams[@]}" --threads 16
    if [[ $? -eq 0 ]]; then
        touch "$output_dir/finish_list/HPF_MC.txt"
        echo "  Done: HPF_MC.bam"
    fi
else
    echo "HPF_MC.bam already exists, skipping."
fi

echo "Merge complete. Output saved to $output_dir."
