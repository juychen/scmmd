#!/bin/bash
allcools merge-allc --chrom_size_path /home/junyichen/code/whole_mouse_brain/wmb/files/mm10.main.chrom.sizes --allc_paths /data2st1/junyi/single_allcs/Vascular/*gz --output_path /data2st1/junyi/merged_allcs/Vascular.tsv.gz
allcools merge-allc --chrom_size_path /home/junyichen/code/whole_mouse_brain/wmb/files/mm10.main.chrom.sizes --allc_paths /data2st1/junyi/single_allcs/OPC-Oligo/*gz --output_path /data2st1/junyi/merged_allcs/OPC-Oligo.tsv.gz
allcools merge-allc --chrom_size_path /home/junyichen/code/whole_mouse_brain/wmb/files/mm10.main.chrom.sizes --allc_paths /data2st1/junyi/single_allcs/Astro-Epen/*gz --output_path /data2st1/junyi/merged_allcs/Astro-Epen.tsv.gz
allcools merge-allc --chrom_size_path /home/junyichen/code/whole_mouse_brain/wmb/files/mm10.main.chrom.sizes --allc_paths /data2st1/junyi/single_allcs/Neuron/*gz --output_path /data2st1/junyi/merged_allcs/Neuron.tsv.gz
allcools merge-allc --chrom_size_path /home/junyichen/code/whole_mouse_brain/wmb/files/mm10.main.chrom.sizes --allc_paths /data2st1/junyi/single_allcs/Immune/*gz --output_path /data2st1/junyi/merged_allcs/Immune.tsv.gz
