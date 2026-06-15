#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate tobias
ulimit -n 65534

cd /data2st2/junyi/output/atac1112/tobiasbam
folder1="/data2st2/junyi/output/atac1112/tobiasbam/MC_PFC/corrected/"
folder2="/data2st2/junyi/output/atac1112/tobiasbam/MW_PFC/corrected/"

folderTF=""

# Define the output folder for merged files
output_dir="/data2st2/junyi/output/atac1112/tobiasbam/PFC_Jaspar26/"
mkdir -p "$output_dir"  # Create output directory if it doesn't exist
# Iterate through all .bam files in the first folder
for bam_file in "$folder1"/*footprints.bw; do
    # Extract the base filename (without path)
    filename=$(basename "$bam_file")
    correctname="${bam_file//_footprints/_corrected}" 
    new_name="${filename//MC/MW}"    # Replace all occurrences of "MW" with "MC"
    correctnameMW="${correctname//MC/MW}"

    for bed_bind in /data2st2/junyi/output/atac1112/tobiasbam/PFC_Jaspar26/$filename/*/beds/*_all.bed; do
        # Check if the same filename exists in the other two folders
        if [[ -f "$folder2/$new_name" ]]; then
            echo "Printagg $filename..."
            echo $folder1/$filename
            echo $folder2/$new_name
            echo $bed_bind
            # Merge the three files into one
            # if the finish file is not in the folder
            ctname=$(basename "${filename%_footprints.bw}")
            TOBIAS PlotAggregate --TFBS $bed_bind  \
            --signals $correctname $correctnameMW \
            --output $bed_bind\.pdf --share_y both --plot_boundaries --signal-on-x \
            --output-txt $bed_bind\_signal.txt
            # if samtools finished successfully
            echo "Finished processing "$folder1/$filename" "$folder2/$new_name" "
            # Index the merged file (optional)
            # samtools index "$output_dir/$filename"
            # If finish the write a text to '/data2st2/junyi/output/atac1112/tobiasbam/MW_PFC/finish_list'
    else
        echo "Skipping $filename (not found in all folders)"
    fi
    done
done

echo "Merging complete. Output saved to $output_dir."