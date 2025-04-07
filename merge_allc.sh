#!/bin/bash
# allcools merge-allc --chrom_size_path /home/junyichen/code/whole_mouse_brain/wmb/files/mm10.main.chrom.sizes --allc_paths /data2st1/junyi/single_allcs/PFC/Vascular/*gz --output_path /data2st1/junyi/merged_allcs/PFC/Vascular.tsv.gz
# allcools merge-allc --chrom_size_path /home/junyichen/code/whole_mouse_brain/wmb/files/mm10.main.chrom.sizes --allc_paths /data2st1/junyi/single_allcs/PFC/OPC-Oligo/*gz --output_path /data2st1/junyi/merged_allcs/PFC/OPC-Oligo.tsv.gz
# allcools merge-allc --chrom_size_path /home/junyichen/code/whole_mouse_brain/wmb/files/mm10.main.chrom.sizes --allc_paths /data2st1/junyi/single_allcs/PFC/Astro-Epen/*gz --output_path /data2st1/junyi/merged_allcs/PFC/Astro-Epen.tsv.gz
# allcools merge-allc --chrom_size_path /home/junyichen/code/whole_mouse_brain/wmb/files/mm10.main.chrom.sizes --allc_paths /data2st1/junyi/single_allcs/PFC/Neuron/*gz --output_path /data2st1/junyi/merged_allcs/PFC/Neuron.tsv.gz
# allcools merge-allc --chrom_size_path /home/junyichen/code/whole_mouse_brain/wmb/files/mm10.main.chrom.sizes --allc_paths /data2st1/junyi/single_allcs/PFC/Immune/*gz --output_path /data2st1/junyi/merged_allcs/PFC/Immune.tsv.gz

# allcools merge-allc --chrom_size_path /home/junyichen/code/whole_mouse_brain/wmb/files/mm10.main.chrom.sizes --allc_paths /data2st1/junyi/single_allcs/AMY/Vascular/*gz --output_path /data2st1/junyi/merged_allcs/AMY/Vascular.tsv.gz
# allcools merge-allc --chrom_size_path /home/junyichen/code/whole_mouse_brain/wmb/files/mm10.main.chrom.sizes --allc_paths /data2st1/junyi/single_allcs/AMY/OPC-Oligo/*gz --output_path /data2st1/junyi/merged_allcs/AMY/OPC-Oligo.tsv.gz
# allcools merge-allc --chrom_size_path /home/junyichen/code/whole_mouse_brain/wmb/files/mm10.main.chrom.sizes --allc_paths /data2st1/junyi/single_allcs/AMY/Astro-Epen/*gz --output_path /data2st1/junyi/merged_allcs/AMY/Astro-Epen.tsv.gz
# allcools merge-allc --chrom_size_path /home/junyichen/code/whole_mouse_brain/wmb/files/mm10.main.chrom.sizes --allc_paths /data2st1/junyi/single_allcs/AMY/Neuron/*gz --output_path /data2st1/junyi/merged_allcs/AMY/Neuron.tsv.gz
# allcools merge-allc --chrom_size_path /home/junyichen/code/whole_mouse_brain/wmb/files/mm10.main.chrom.sizes --allc_paths /data2st1/junyi/single_allcs/AMY/Immune/*gz --output_path /data2st1/junyi/merged_allcs/AMY/Immune.tsv.gz

# allcools merge-allc --chrom_size_path /home/junyichen/code/whole_mouse_brain/wmb/files/mm10.main.chrom.sizes --allc_paths /data2st1/junyi/single_allcs/HIP/Vascular/*gz --output_path /data2st1/junyi/merged_allcs/HIP/Vascular.tsv.gz
# allcools merge-allc --chrom_size_path /home/junyichen/code/whole_mouse_brain/wmb/files/mm10.main.chrom.sizes --allc_paths /data2st1/junyi/single_allcs/HIP/OPC-Oligo/*gz --output_path /data2st1/junyi/merged_allcs/HIP/OPC-Oligo.tsv.gz
# allcools merge-allc --chrom_size_path /home/junyichen/code/whole_mouse_brain/wmb/files/mm10.main.chrom.sizes --allc_paths /data2st1/junyi/single_allcs/HIP/Astro-Epen/*gz --output_path /data2st1/junyi/merged_allcs/HIP/Astro-Epen.tsv.gz
source /home/junyichen/anaconda3/etc/profile.d/conda.sh

conda activate allcools
cd /data2st1/junyi/single_allcs/HIP/Neuron
allcools merge-allc --chrom_size_path /home/junyichen/code/whole_mouse_brain/wmb/files/mm10.main.chrom.sizes --allc_paths /data2st1/junyi/single_allcs/HIP/Neuron/hip_neuron.txt --output_path /data2st1/junyi/merged_allcs/HIP/Neuron.tsv.gz
# allcools merge-allc --chrom_size_path /home/junyichen/code/whole_mouse_brain/wmb/files/mm10.main.chrom.sizes --allc_paths /data2st1/junyi/single_allcs/HIP/Immune/*gz --output_path /data2st1/junyi/merged_allcs/HIP/Immune.tsv.gz

