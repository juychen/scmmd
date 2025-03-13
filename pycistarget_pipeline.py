# %%
import pandas as pd
import re
import os


# %%
filename = "/data2st1/junyi/output/motif/AMY_Neuron_MC_wilcoxon.csv"

# %%
df_dar  = pd.read_csv(filename,index_col=0)

# %%
experimentname = re.split(r'[./]', filename)[-2]
experimentname
folder_name = os.path.dirname(filename)

# %%
df_dar_filtered = df_dar[(df_dar['pvals']<0.05) & (df_dar['logfoldchanges']>0) ]

# %%
df_dar_filtered.sort_values(by='logfoldchanges',ascending=False,inplace=True)

# %%
# For liftover
df_dar_filtered.names.str.split(r'[ ,!\-;:|]',expand=True).to_csv(f"{folder_name}/{experimentname}_DAR.bed",header=False,index=False,sep="\t")

# %%
import pyranges as pr
import os
path_to_region_sets = '/data2st1/junyi/output/motif'
region_sets_files = ['AMY_Neuron_MC_wilcoxon_DAR_lift.bed','AMY_Neuron_MC_wilcoxon_DAR_lift.bed']
region_sets = {x.replace('.bed', ''):pr.read_bed(os.path.join(path_to_region_sets, x)) for x in region_sets_files}


# %%
from pycistarget.motif_enrichment_cistarget import *


# %%
region_sets

# %%
# !pycistarget cistarget --cistarget_db_fname '/data2st1/junyi/scenic/mouse/motif/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather' --bed_fname '/data2st1/junyi/output/motif/AMY_Neuron_MC_wilcoxon_DAR_lift.bed' \
# --species 'mus_musculus' \
# --auc_threshold 0.005 \
# --nes_threshold 3.0 \
# --rank_threshold 0.05 \
# --path_to_motif_annotations '/data2st1/junyi/scenic/mouse/motif/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl' \
# --output_folder '/data2st1/junyi/output/motif/' \
# --write_html


# %%
# cistarget_dict = run_cistarget_command(ctx_db = '/data2st1/junyi/scenic/mouse/motif/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather',
#                                                       region_sets = region_sets,
#                                                       specie = 'mus_musculus',
#                                                       auc_threshold = 0.005,
#                                                       nes_threshold = 3.0,
#                                                       rank_threshold = 0.05,
#                                                       annotation = ['Direct_annot', 'Orthology_annot'],
#                                                       annotation_version = 'v10nr_clust',
#                                                       path_to_motif_annotations = '/data2st1/junyi/scenic/mouse/motif/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl',
#                                                       n_cpu = 4,
#                                                       _temp_dir='/data2st1/junyi/output/motif')

regions: pr.PyRanges = pr.read_bed('/data2st1/junyi/output/motif/'+region_sets_files[0], as_df=False)

ctx_db = cisTargetDatabase(
    fname='/data2st1/junyi/scenic/mouse/motif/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather',
    region_sets=regions,
)




