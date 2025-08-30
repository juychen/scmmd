#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate tobias
ulimit -n 65534

cd /data2st2/junyi/output/atac0627/tobiasbam
folder1="/data2st2/junyi/output/atac0627/tobiasbam/MC_PFC/corrected/"
folder2="/data2st2/junyi/output/atac0627/tobiasbam/MW_PFC/corrected/"

# Define the output folder for merged files
output_dir="/data2st2/junyi/output/atac0627/tobiasbam/PFC"
mkdir -p "$output_dir"  # Create output directory if it doesn't exist

# Iterate through all .bam files in the first folder
for bam_file in "$folder1"/*footprints.bw; do
    # Extract the base filename (without path)
    filename=$(basename "$bam_file")
    new_name="${filename//MC/MW}"    # Replace all occurrences of "MW" with "MC"

    # Check if the same filename exists in the other two folders
    if [[ -f "$folder2/$new_name" ]]; then
        echo "Merging $filename..."
        echo $folder1/$filename
        echo $folder2/$new_name
        # Merge the three files into one
        # if the finish file is not in the folder
        # 
        mkdir -p "$output_dir/finish_list"
        if [[ ! -f "$output_dir/finish_list/$filename.txt" ]]; then
            mkdir -p "$output_dir/$filename"
            TOBIAS BINDetect --motifs /data2st1/junyi/scenic/mouse/motif/merged_cluster/direct_im_motifs.jaspar \
            --signals "$folder1/$filename" "$folder2/$new_name" \
            --outdir "$output_dir/$filename" --cond_names MC MW --cores 32 \
            --genome /data2st1/junyi/ref/GRCm38.p6.genome.fa \
            --peaks /data2st1/junyi/output/atac0627/cCRE/peak.bed \
            # if samtools finished successfully
            if [[ $? -eq 0 ]]; then
                touch "$output_dir/finish_list/$filename.txt"
            fi
            echo "Finished processing "$folder1/$filename" "$folder2/$new_name" "
        fi
        # Index the merged file (optional)
        # samtools index "$output_dir/$filename"
        # If finish the write a text to '/data2st2/junyi/output/atac0627/tobiasbam/MW_PFC/finish_list'
    else
        echo "Skipping $filename (not found in all folders)"
    fi
done

echo "Merging complete. Output saved to $output_dir."
