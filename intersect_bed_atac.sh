#!/bin/bash
source activate allcools
bedtools intersect -a /data8/junyi/generegion_vM33/region_anotation_typeonly.bed -b /data8/hannan/snATAC_process/snATAC_01_bam/WT_W26/outs/fragments.tsv.gz -c > /data8/junyi/pseudobulk_methlylation/atac_cre_counts.bed
