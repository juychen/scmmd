#!/bin/bash

# Check if the correct number of arguments is provided
# if [ "$#" -ne 2 ]; then
#     echo "Usage: $0 <input_folder> <output_folder>"
#     exit 1
# fi

# Assign input and output folders
input_folder="/data2st1/junyi/atac_bigbed/"
output_folder="/data2st1/junyi/atac_bigbed/peak/"

conda activate allcools
# Check if the input folder exists
if [ ! -d "$input_folder" ]; then
    echo "Error: Input folder '$input_folder' does not exist."
    exit 1
fi

# Create the output folder if it doesn't exist
mkdir -p "$output_folder"

# Check if bigWigToBedGraph is installed
if ! command -v bigWigToBedGraph &> /dev/null; then
    echo "Error: 'bigWigToBedGraph' tool not found. Please install UCSC tools."
    exit 1
fi

# Loop through all .bw files in the input folder
for bw_file in "$input_folder"/*.bedGraph; do
    # Get the base name of the file (without extension)
    base_name=$(basename "$bw_file" .bedGraph)

    # Define the output BedGraph file path
    bedgraph_file="$output_folder/$base_name.bed"

    # Convert BigWig to BedGraph
    echo "Converting $bw_file to $bedgraph_file..."
    macs2 bdgpeakcall -i "$bw_file" -o "$bedgraph_file" --cutoff 5

    # Check if the conversion was successful
    if [ $? -eq 0 ]; then
        echo "Saved BedGraph file: $bedgraph_file"
    else
        echo "Error converting $bw_file"
    fi
done

echo "Conversion complete."