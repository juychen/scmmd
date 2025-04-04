# %%
from ALLCools.mcds import RegionDS
from ALLCools.dmr import call_dms, call_dmr
import glob as glob
import pandas as pd

# %%
cell_types = ['Neuron','Astro-Epen','Vascular','Immune','OPC-Oligo']

# %%
mc_bulk_dir = '/data2st1/junyi/merged_allcs/PFC/'
allc_paths = glob.glob(f'{mc_bulk_dir}/*.tsv.gz')

# %%
celltype = cell_types[1]

# %%
# 3rd column is the orginal id in the wilcox table
df_dar = pd.read_csv(f'/data2st1/junyi/output/motif/ALL_{celltype}_wilcoxon_DAR_lifted_sorted_gene.bed',sep='\t',header=None,index_col=3)

# %%
df_dar.columns = ['chr','start','end','chrmm10','startmm10','endmm10','score','strand','gene_name','gene_id','region','distance']

# %%
df_dar

# %%
df_pval = pd.read_csv(f'/data2st1/junyi/output/motif/ALL_{celltype}_wilcoxon.csv',index_col=0)

# %%
df_concatted= pd.concat([df_pval.loc[df_dar.index],df_dar],axis=1)

# %%
df_dar_selected = df_concatted[(df_concatted.pvals<0.05) & (df_concatted.logfoldchanges>0.1)& (df_concatted.pct_nz_group>0.1)]

# %%
df_dar_selected['regionmm10'] = df_dar_selected['chr'] + ':' + df_dar_selected['start'].astype(str) + '-' + df_dar_selected['end'].astype(str)

# %%
df_dar_selected['region_extended'] = df_dar_selected['chr'] + ':' + (df_dar_selected['start']-5000).astype(str) + '-' + (df_dar_selected['end']+5000).astype(str)

# %%


# %%
df_dar_selected.sort_values(by='scores',ascending=False).head(10)

# %%
allc_table= {}

for allc_path in cell_types:
    allc_table["/data2st1/junyi/merged_allcs/PFC/"+allc_path] = mc_bulk_dir+allc_path+'.tsv.gz'

# %%
allc_table

# %%
# make a dict, key is sample name, value is allc path
# import pathlib
# allc_table = {
#     allc_path.name.split('.')[0]: str(allc_path)
#     for allc_path in pathlib.Path(mc_bulk_dir).glob(
#         f'{mc_bulk_dir}/*.tsv.gz')
# }
samples = list(allc_table.keys())
allc_paths = list(allc_table.values())

chrom_size_path = '/home/junyichen/code/whole_mouse_brain/wmb/files/mm10.main.chrom.sizes'
output_dir = '/data2st1/junyi/output/dmr/test_PFC'


# %%
input_region = df_dar_selected.sort_values(by='scores',ascending=False).head(10)
input_region.sort_values(['chr','start','end'],inplace=True)
input_region.drop_duplicates(subset=['chr','start','end'],keep='first',inplace=True)
input_region = input_region.regionmm10.values.tolist()  

# %%
call_dms(
    output_dir=output_dir,
    allc_paths=allc_paths,
    samples=samples,
    chrom_size_path=chrom_size_path,
    cpu=45,
    max_row_count=50,
    n_permute=3000,
    min_pvalue=0.01,
    # here we just calculate some small regions for demo
    # do not provide region parameter if you want to run DMR calling for the whole genome
    # This parameter can also be used for call DMR/DMS in specific region of interest
    region=['chr7:4106641-4137142'])


# %%
call_dmr(output_dir=output_dir,
         p_value_cutoff=0.001,
         frac_delta_cutoff=0.3,
         max_dist=250,
         residual_quantile=0.7,
         corr_cutoff=0.3,
         cpu=30)

# %%
dmr_ds = RegionDS.open(output_dir)
dmr_ds


# %%
dmr_ds['dmr_da_frac'].to_numpy()[3]

# %%
dmr_ds['dmr_da_frac'].to_numpy()[1]

# %%
dmr_ds['dmr_da_frac'].to_numpy()

# %%
dmr_ds['sample'].to_numpy()

# %%
dmr_ds['dmr_chrom'].to_numpy()

# %%
dmr_ds['dmr_start'].to_numpy()

# %%
dmr_ds['dmr_end'].to_numpy()

# %%
dmr_ds['dmr_state'].to_numpy()[1]

# %%
dmr_ds['dmr_state'].to_numpy()[3]

# %%
dmr_ds['dmr_state'].to_numpy().sum(axis=1)

# %%
dmr_ds['dmr_da_frac'].to_numpy()

# %%
#dmr_ds = RegionDS.open(output_dir, region_dim='dms')

# %%
dmr_ds['dms']


