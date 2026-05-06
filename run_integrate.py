# %%
import matplotlib.pyplot as plt
import anndata as adata
import scanpy as sc
from wmb import cemba, mm10
from ALLCools.mcds import MCDS
from ALLCools.clustering import significant_pc_test
from ALLCools.plot import *
import pandas as pd
from harmonypy import run_harmony

# %%
adata_sc = adata.read_h5ad("/data1st1/hydata/pfc_3tech_int_byID.h5ad")

# %%
adata_sc.var.columns

# %%
adata_sc.obs.columns

# %%
adata_sc.obsm.keys()

# %%
hue = 'total_counts'
fig, axes = plot_decomp_scatters(adata_sc,
                                 n_components=10,
                                 hue=hue,
                                 hue_quantile=(0.25, 0.75),
                                 nrows=5,
                                 ncols=5)


# %%
fig, axes = plt.subplots(figsize=(8, 4), dpi=300, ncols=2)
ax = axes[0]
categorical_scatter(ax=ax, data=adata_sc, hue='louvain', palette='tab20')
ax = axes[1]
categorical_scatter(ax=ax, data=adata_sc, hue='pred_mwb', palette='tab20')


# %% [markdown]
# # Prepare L4 region level methlylation data and load the gene 2k flaking bases

# %%
var_dim = 'geneslop2k-vm23'
mcds = MCDS.open('/data2st1/junyi/methlyatlas/mCseq/data.nemoarchive.org/biccn/grant/u19_cemba/ecker/epigenome/cellgroup/mCseq3/mouse/processed/counts/CEMBA.snmC.L4RegionAgg.zarr',
                   obs_dim='L4Region',var_dim=var_dim)


# %%
mcds

# %%
gene_meta_path = '/home/junyichen/code/whole_mouse_brain/wmb/files/modified_gencode.vM23.primary_assembly.annotation.gene.flat.tsv.gz'


# %%
gene_meta = pd.read_csv(gene_meta_path, index_col='gene_id', sep='\t')


# %%
gene_meta.head()

# %%
import pybedtools

genes_to_skip = set()
chrom_to_remove = ['chrM']


# skip smaller genes mostly covered by a larger gene, e.g., a miRNA within a protein coding gene.
# F=0.9 means > 90% of gene_b is overlapped with gene_a, in this case, we only keep gene_a for DMG test
gene_bed = pybedtools.BedTool.from_dataframe(
    gene_meta.reset_index()[['chrom', 'start', 'end', 'gene_id']])
mapped_bam = gene_bed.map(b=gene_bed, c=4, o='distinct', F=0.9)
for _, (*_, gene_a, gene_b_str) in mapped_bam.to_dataframe().iterrows():
    for gene_b in gene_b_str.split(','):
        if gene_b != gene_a:
            genes_to_skip.add(gene_b)

# remove certain chromosomes
genes_to_skip |= set(gene_meta.index[gene_meta['chrom'].isin(chrom_to_remove)])
use_features = gene_meta.index[~gene_meta.index.isin(genes_to_skip)]
print(f'{use_features.size} features remained')


# %%
min_cov = 5
mcds.add_feature_cov_mean()

feature_cov_mean = mcds.coords[f'{var_dim}_cov_mean'].to_pandas()
use_features &= feature_cov_mean[feature_cov_mean > min_cov].index

print(f'{use_features.size} features remained')


# %%
use_features

# %%
mcds.filter_feature_by_cov_mean(min_cov=min_cov)


# %%
mcds

# %%
mcds.add_mc_frac(normalize_per_cell=True, clip_norm_value=10)


# %%
mcds

# %%
mcds = mcds[['geneslop2k-vm23_da_frac']]
mcds['geneslop2k-vm23_da_frac'] = mcds['geneslop2k-vm23_da_frac'].astype('float32')


# %%
mcds

# %%
mcds

# %%
mcds.write_dataset('output/geneslop2k-vm23_frac.mcds', var_dims=['geneslop2k-vm23'],mode='w')


# %%
use_gene_meta = gene_meta.loc[use_features]
use_gene_meta.to_csv('output/GeneMetadata.csv.gz')


# %%
mcg = mcds.get_adata(mc_type='CGN',var_dim='geneslop2k-vm23',obs_dim='L4Region')

# %%
mcg


# %%
mcg.obs

# %%
mch = mcds.get_adata(mc_type='CHN',var_dim='geneslop2k-vm23',obs_dim='L4Region')

# %%
mch

# %% [markdown]
# # Prepare L4 region level metadata

# %%
meta_data = pd.read_csv('/data2st1/junyi/methlyatlas/mCseq/CEMBA.mC.Metadata/CEMBA.mC.Metadata.csv')

# %%
meta_data.head()

# %%
meta_data.columns

# %% [markdown]
# # Group metat for cell group level
# - some are median 
# - some are sum

# %%
result = meta_data.groupby('CellGroup').agg(lambda x: x.value_counts().idxmax())[['CEMBARegion', 'MajorRegion', 'SubRegion', 'SubClass',
       'Class', 'NeuroTransmitters']]


# %%
result

# %%
meta_median = meta_data.groupby('CellGroup').median()[[ 'mCCCFrac', 'mCGFrac', 'mCHFrac','PlateNormCov','FinalmCReads','InputReads']]
meta_median


# %%
meta_median.columns = [col + '_median' for col in meta_median.columns]

# %%
meta_median

# %%
meta_sum = meta_data.groupby('CellGroup').sum()[['FinalmCReads','InputReads']]
meta_sum.columns = [col + '_sum' for col in meta_sum.columns]
meta_sum

# %%
meta_group = pd.read_csv('/data2st1/junyi/methlyatlas/mCseq/CEMBA.mC.CellGroup.Coordinates/CEMBA.mC.CellGroup.Coordinates.csv.gz',index_col=0)

# %%
meta_group

# %%
if result.index.equals(meta_median.index) and result.index.equals(meta_sum.index) and result.index.equals(meta_group.index):
    meta_cell_grop = pd.concat([result,meta_median,meta_sum,meta_group],axis=1)
else:
    print('index not match')

# %%
meta_cell_grop.to_csv('output/meta_cell_group.csv.gz')

# %%
mcg.obs = mcg.obs.join(meta_cell_grop)

# %%
mcg.obsm['X_umap'] = mcg.obs[['mc_all_umap_0', 'mc_all_umap_1']].values


# %%
mch.obs = mch.obs.join(meta_cell_grop)

# %%
mch.obsm['X_umap'] = mch.obs[['mc_all_umap_0', 'mc_all_umap_1']].values


# %%
sc.pl.umap(mcg, color=['Class'], ncols=2)

# %%
sc.pl.umap(mch, color=['Class'], ncols=2)

# %%


# %% [markdown]
# # Prepare the mcds file will 100k bin 

# %%
var_dim = 'chrom100k'
mcds_bin = MCDS.open('/data2st1/junyi/methlyatlas/mCseq/data.nemoarchive.org/biccn/grant/u19_cemba/ecker/epigenome/cellgroup/mCseq3/mouse/processed/counts/CEMBA.snmC.L4RegionAgg.zarr',
                   obs_dim='L4Region',var_dim='chrom100k',
                    use_obs=meta_cell_grop.index
                   )

total_feature = mcds_bin.get_index(var_dim).size
mcds_bin


# %%
mcds_bin.add_cell_metadata(meta_cell_grop)


# %%
mcds_bin.add_feature_cov_mean(var_dim=var_dim)


# %%

# feature cov cutoffs
min_cov = 500
max_cov = 3000
black_list_fraction = 0.2
black_list_path = mm10.ENCODE_BLACKLIST_PATH
exclude_chromosome = ['chrM']
# HVF
mch_pattern = 'CHN'
mcg_pattern = 'CGN'
n_top_feature = 20000

# PC cutoff
pc_cutoff = 0.1

# KNN
knn = -1  # -1 means auto determine

# Leiden
resolution = 1


# filter by coverage - based on the distribution above
mcds_bin = mcds_bin.filter_feature_by_cov_mean(
    min_cov=min_cov,  # minimum coverage
    max_cov=max_cov  # maximum coverage
)

# remove blacklist regions
mcds_bin = mcds_bin.remove_black_list_region(
    black_list_path=black_list_path,
    f=black_list_fraction  # Features having overlap > f with any black list region will be removed.
)

# remove chromosomes
mcds_bin = mcds_bin.remove_chromosome(exclude_chromosome)

# %%
mcds_bin

# %%
load = True
mcds_bin.add_mc_frac(
normalize_per_cell=True,  # after calculating mC frac, per cell normalize the matrix
    clip_norm_value=10  # clip outlier values above 10 to 10
)
obs_dim = 'L4Region'

# load only the mC fraction matrix into memory so following steps is faster
# Only load into memory when your memory size is enough to handle your dataset
if load and (mcds_bin.get_index(obs_dim).size < 20000):
    mcds_bin[f'{var_dim}_da_frac'].load()


# %%
mch_bin_hvf = mcds_bin.calculate_hvf_svr(var_dim=var_dim,
                                 mc_type=mcg_pattern,
                                 n_top_feature=n_top_feature,
                                 plot=True)

# %%
mcg_bin_adata = mcds_bin.get_adata(mc_type=mcg_pattern,
                           var_dim=var_dim,
                           select_hvf=True)
mcg_bin_adata


# %%
mcg_bin_adata.obs.head()

# %%
from ALLCools.clustering import tsne, significant_pc_test, log_scale
from ALLCools.plot import *

# %%
mcg_bin_adata

# %%
log_scale(mcg_bin_adata)


# %%
mcg_col_name = 'mCGFrac_median'  # Name may change

sc.tl.pca(mcg_bin_adata)
cg_n_components = significant_pc_test(mcg_bin_adata)
fig, axes = plot_decomp_scatters(mcg_bin_adata,
                                 n_components=cg_n_components,
                                 hue=mcg_col_name,
                                 hue_quantile=(0.25, 0.75),
                                 nrows=3,
                                 ncols=5)

# %%
mch_hvf = mcds_bin.calculate_hvf_svr(var_dim=var_dim,
                                 mc_type=mch_pattern,
                                 n_top_feature=n_top_feature,
                                 plot=True)


# %%
mch_bin_adata = mcds_bin.get_adata(mc_type=mch_pattern,
                           var_dim=var_dim,
                           select_hvf=True)
mch_bin_adata


# %%
mch_bin_adata

# %%
mch_col_name = "mCHFrac_median"
sc.tl.pca(mch_bin_adata)
ch_n_components = significant_pc_test(mch_bin_adata)
fig, axes = plot_decomp_scatters(mch_bin_adata,
                                 n_components=ch_n_components,
                                 hue=mch_col_name,
                                 hue_quantile=(0.25, 0.75),
                                 nrows=3,
                                 ncols=5)

# %%
import numpy as np
ch_pcs = mch_bin_adata.obsm['X_pca'][:, :ch_n_components]
cg_pcs = mcg_bin_adata.obsm['X_pca'][:, :cg_n_components]

cg_pcs = cg_pcs / cg_pcs.std()
ch_pcs = ch_pcs / ch_pcs.std()

# total_pcs
total_pcs = np.hstack([ch_pcs, cg_pcs])

# make a copy of adata, add new pcs
# this is suboptimal, will change this when adata can combine layer and X in the future
adata_merged = mch_bin_adata.copy()
adata_merged.obsm['X_pca'] = total_pcs


# %%
adata_merged

# %%
if knn == -1:
    knn = max(15, int(np.log2(adata_merged.shape[0])*2))


# %%
sc.pp.neighbors(adata_merged, n_neighbors=knn)


# %%
sc.tl.leiden(adata_merged, resolution=resolution)


# %%
sc.tl.umap(adata_merged)



# %%
fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
_ = categorical_scatter(data=adata_merged,
                        ax=ax,
                        coord_base='umap',
                        hue='leiden',
                        text_anno='leiden',
                        show_legend=True)


# %%
fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
_ = categorical_scatter(data=adata_merged,
                        ax=ax,
                        coord_base='umap',
                        hue='Class',
                        text_anno='Class',
                        show_legend=True)


# %%
fig, ax = plt.subplots(figsize=(4, 4), dpi=300)
_ = categorical_scatter(data=adata_merged,
                        ax=ax,
                        coord_base='umap',
                        hue='leiden',
                        text_anno='leiden',
                        show_legend=True)

# %%
pd.DataFrame(set(adata_merged.obs['SubClass'])).to_csv('output/subclass.csv.gz')

# %%
interactive_scatter(data=mcg,
                    hue='Class',
                    coord_base='umap',
                    max_points=3000)

# %%
interactive_scatter(data=adata_merged,
                    hue='Class',
                    coord_base='umap',
                    max_points=3000)

# %%
adata_sc.obs.columns

# %%
map_subclass_class = dict(zip(adata_merged.obs['SubClass'], adata_merged.obs['Class']))

# %%
cell_class = adata_sc.obs['pred_mwb'].str.replace(r'^\d+\s', '')
convert_class = []
for i in range(0,len(cell_class)):
    cc = cell_class[i] 
    if cc not in map_subclass_class:
        convert_class.append(cc)
    else:
        convert_class.append(map_subclass_class[cc])   


# %%
adata_sc.obs['pred_Class'] = convert_class

# %%
interactive_scatter(data=adata_sc,
                    hue='pred_mwb',
                    coord_base='umap',
                    max_points=3000)

# %%
interactive_scatter(data=adata_sc,
                    hue='pred_Class',
                    coord_base='umap',
                    max_points=3000)

# %%
interactive_scatter(data=adata_merged,
                    hue='leiden',
                    coord_base='umap',
                    max_points=3000)

# %%
sc.pl.umap(mcg, color=['Class'], ncols=2)

# %%
adata_merged.write_h5ad('output/Brain.chrom100k-clustering.h5ad')
adata_merged


# %%
adata_merged.obs.to_csv('output/Brain.ClusteringResults.csv.gz')
adata_merged.obs.head()


# %%
mcds

# %%
sc.pp.highly_variable_genes(mcg, n_top_genes=2000, n_bins=20)

# %%
sc.pp.highly_variable_genes(mch, n_top_genes=2000, n_bins=20)

# %%
adata_raw = sc.read_10x_mtx("/data1st1/seq_from_our_lab/NovaSeq_mmBrain_X101SC24104268-Z01/X101SC24104268-Z01-J001/seeksoul_output/W1/step3/filtered_feature_bc_matrix")

# %%
adata_raw.var

# %%
n_genes =5000

# %%
adata_sc.var['gane_name'] = adata_sc.var.index

# %%
adata_sc.var

# %%
var_merged = adata_raw.var.merge(adata_sc.var,how='outer',left_index=True,right_index=True,indicator=True)
var_merged

# %%
adata_sc.var.index

# %%
adata_sc.var = var_merged.loc[adata_sc.var.index]

# %%
adata_sc.var.set_index('gene_ids',inplace=True)

# %%
adata_atac_gene = adata.read_h5ad('/data2st1/junyi/snATAC_02_matrix/WT_PFC_W26_geneactivity_n.h5ad')
adata_atac_bin = adata.read_h5ad('/data2st1/junyi/snATAC_02_matrix/WT_PFC_W26_5kb.h5ad')

# %%
pcs = pd.read_csv('/data2st1/junyi/snATAC_02_matrix/combined_atac_QC_reduction_pcs.csv',index_col=0)

# %%
pcs.index

# %%
len(set("WT_"+adata_atac_bin.obs.index).intersection(set(pcs.index)))

# %%
adata_atac_bin.obs.index = "WT_"+adata_atac_bin.obs.index

# %%
adata_atac_gene.obs.index = "WT_"+adata_atac_gene.obs.index

# %%
adata_atac_gene

# %%
merge_barcode = adata_atac_gene.obs.merge(pcs,left_index=True,right_index=True,how='inner')

# %%
adata_atac_selected = adata_atac_gene[merge_barcode.index]

# %%
adata_atac_selected.obsm['X_pca'] = pcs.loc[merge_barcode.index].iloc[:,1:].values
# Ignore the first pc 

# %%
adata_atac_selected.obsm['X_pca'].shape

# %%
adata_atac_selected

# %%
adata_atac_selected.obs.head()

# %%
df_atac_umapdf_atac_umap = pd.read_csv('/data2st1/junyi/snATAC_02_matrix/combined_atac_QC_reduction_umap.csv',index_col=0)
adata_atac_selected.obsm['X_umap_seurat'] = df_atac_umapdf_atac_umap.loc[merge_barcode.index].values

# %%
def dump_embedding(adata, name, n_dim=2):
    # put manifold coordinates into adata.obs
    for i in range(n_dim):
        adata.obs[f'{name}_{i}'] = adata.obsm[f'X_{name}'][:, i]
    return


# %%
sc.pp.neighbors(adata_atac_selected)
sc.tl.leiden(adata_atac_selected, resolution=1)
sc.tl.umap(adata_atac_selected)
dump_embedding(adata_atac_selected, 'umap')


# %%
dump_embedding(adata_atac_selected, 'umap')
dump_embedding(adata_atac_selected, 'umap_seurat')


# %%
fig, axes = plt.subplots(figsize=(12, 4), dpi=250, ncols=2)

ax = axes[0]
_ = categorical_scatter(data=adata_atac_selected.obs,
                        ax=ax,
                        coord_base='umap',
                        hue='leiden',
                        palette='tab20',
                        text_anno='leiden')
ax = axes[1]
_ = categorical_scatter(data=adata_atac_selected.obs,
                        ax=ax,
                        coord_base='umap_seurat',
                        hue='leiden',
                        palette='tab20',
                        text_anno='leiden')



# %%
n_components = significant_pc_test(adata_atac_selected)


# %%
fig, axes = plot_decomp_scatters(adata_atac_selected,
                                 n_components=n_components,
                                 hue=None,
                                 hue_quantile=(0.25, 0.75),
                                 nrows=5,
                                 ncols=5)

# %%
n_genes = 5000
sc.pp.highly_variable_genes(adata_atac_selected, n_top_genes=n_genes, n_bins=20) 

# %%
adata_atac_selected.var['gene_id'] = adata_atac_selected.var.index.str.split('.').str[0]

# %%
adata_atac_selected.var.set_index('gene_id',inplace=True)

# %%
sc_hvf_list = adata_sc.var['variances_norm'].sort_values(ascending=False).dropna()[:n_genes].index


# %%
mcg_hvf_list = mcg.var['dispersions_norm'].sort_values(ascending=False).dropna()[:n_genes].index


# %%
mch_hvf_list = mch.var['dispersions_norm'].sort_values(ascending=False).dropna()[:n_genes].index

# %%
adata_atac_selected.var['dispersions_norm'].sort_values(ascending=False).dropna()[:n_genes].index

# %%
atac_hvf_list = adata_atac_selected.var['dispersions_norm'].sort_values(ascending=False).dropna()[:n_genes].index


# %%
hvfs = mcg_hvf_list.intersection(sc_hvf_list).intersection(mcg_hvf_list).intersection(atac_hvf_list)

# %%
len(hvfs)

# %%
sc_hvf_adata = adata_sc[:,hvfs].copy()
sc.pp.scale(sc_hvf_adata)
sc_hvf_adata

# %%
mcg_hvf_adata = mcg[:,hvfs].copy()
sc.pp.scale(mcg_hvf_adata)
mcg_hvf_adata.X *=-1

# %%
mch_hvf_adata = mch[:,hvfs].copy()
sc.pp.scale(mch_hvf_adata)
mch_hvf_adata.X *=-1

# %%
atac_hvf_adata = adata_atac_selected[:,hvfs].copy()
sc.pp.scale(atac_hvf_adata)

# %%
adata_concat = sc_hvf_adata.concatenate([
    mcg_hvf_adata,
    mch_hvf_adata,
    atac_hvf_adata
    
],
                                index_unique='-')
adata_concat.obs['batch'] = adata_concat.obs['batch'].map({
    '0': 'sc',
    '1': 'mcg',
    '2': 'mch',
    '3': 'atac'
    })
adata_concat

# %%
sc.pp.pca(adata_concat, n_comps=50)
sc.pl.pca_variance_ratio(adata_concat)
n_components = significant_pc_test(adata_concat, p_cutoff=0.2)


# %%
adata_atac_selected.obs.columns

# %%
adata_concat.obs 

# %%
pd.DataFrame(adata_concat.obs['batch'])

# %%
ho = run_harmony(adata_concat.obsm['X_pca'],
                 meta_data=pd.DataFrame(adata_concat.obs['batch']),
                 vars_use='batch',
                 random_state=0,
                 nclust=100,
                 max_iter_harmony=20)
#adata_concat.obsm['X_pca'] = ho.Z_corr.T


# %%
adata_concat.write_h5ad('output/merged-harmony-clustering.h5ad')

# %%
adata_concat.obsm['X_pca'] = ho.Z_corr.T


# %%
sc.pp.neighbors(adata_concat, n_neighbors=25)
sc.tl.leiden(adata_concat, resolution=1.5)
sc.tl.umap(adata_concat)


# %%
_ = plot_decomp_scatters(adata_concat,
                         n_components=n_components,
                         hue='batch',
                         palette='tab10')

# %%
adata_concat.obs['umap_0'] = adata_concat.obsm['X_umap'][:, 0]
adata_concat.obs['umap_1'] = adata_concat.obsm['X_umap'][:, 1]


# %%
fig, axes = plt.subplots(figsize=(8, 4), dpi=250, ncols=2)

ax = axes[0]
categorical_scatter(ax=ax,
                    data=adata_concat,
                    hue='batch',
                    show_legend=True,
                    max_points=None,
                    s=1)

ax = axes[1]
categorical_scatter(ax=ax,
                    data=adata_concat,
                    hue='pred_Class',
                    show_legend=True,
                    max_points=None,
                    s=1)

# %%
adata_concat.obs.columns

# %%
adata_concat.obs[adata_concat.obs['batch'] == 'mcg'].copy().columns

# %%
fig, axes = plt.subplots(figsize=(8, 4), dpi=250, ncols=2)

data = adata_concat.obs[adata_concat.obs['batch'] == 'mcg'].copy()
#data['MajorType'] = data['SubType'].str.split('_').str[-1]
data2 = adata_concat.obs[adata_concat.obs['batch'] == 'mch'].copy()

ax = axes[0]
categorical_scatter(ax=ax,
                    data=data,
                    show_legend=False,
                    max_points=None,
                    hue='Class',
                    text_anno='Class',
                    s=1)

ax = axes[1]
categorical_scatter(ax=ax,
                    data=adata_concat,
                    hue='pred_Class',
                    show_legend=True,
                    text_anno='pred_Class',

                    max_points=None,
                    s=1)

# %%
adata_sc.obs.columns

# %%
df_pred_class=pd.DataFrame(pd.crosstab(adata_concat.obs['pred_Class'], adata_concat.obs['leiden']).idxmax())

# %%
df_class=pd.DataFrame(pd.crosstab(adata_concat.obs['Class'], adata_concat.obs['leiden']).idxmax())

# %%
df_type_merged = df_pred_class.merge(df_class,how='outer',left_index=True,right_index=True).fillna(method='ffill').fillna(method='bfill')
df_type_merged.columns = ['pred_Class','Class']

# %%
df_type_merged

# %%
lei_map_class = dict(zip(df_type_merged.index,df_type_merged['Class']))
lei_map_pred_class = dict(zip(df_type_merged.index,df_type_merged['pred_Class']))

# %%
adata_concat.obs

# %%
adata_concat.obs

# %%
mapped_class = []
mappled_preclass = []
for inde,row in (adata_concat.obs).iterrows():
    if row['leiden'] not in lei_map_class:
        mapped_class.append(row['leiden'])
        mappled_preclass.append(row['leiden'])
    else:
        mapped_class.append(lei_map_class[row['leiden']])
        mappled_preclass.append(lei_map_pred_class[row['leiden']])  


# %%
adata_concat.obs['mapped_class'] = mapped_class
adata_concat.obs['mappled_preclass'] = mappled_preclass

# %%
interactive_scatter(data=adata_concat,
                    hue='mappled_preclass',
                    coord_base='umap',
                    max_points=20000)

# %%
interactive_scatter(data=adata_concat,
                    hue='mappled_preclass',
                    coord_base='umap',
                    max_points=10000)

# %%
adata_concat.obs

# %%
sc._utils.sanitize_anndata(adata_concat)

# %%
adata_concat.obs['predicted_doublet'].fillna('False',inplace=True)
adata_concat.obs["predicted_doublet"] = adata_concat.obs["predicted_doublet"].astype(bool)


# %%
adata_concat.var['mt-0'] = adata_concat.var["mt-0"].astype(bool)
adata_concat.var['highly_variable-0-0'] = adata_concat.var["highly_variable-0-0"].astype(bool)
adata_concat.var['highly_variable-1-0'] = adata_concat.var["highly_variable-1-0"].astype(bool)
adata_concat.var['highly_variable-0'] = adata_concat.var["highly_variable-0"].astype(bool)


# %%
adata_concat.write_h5ad('output/merged-clustering.h5ad')

# %%
adata_atac_gene = adata.read_h5ad('/data2st1/junyi/snATAC_02_matrix/WT_PFC_W26_geneactivity.h5ad')
adata_atac_bin = adata.read_h5ad('/data2st1/junyi/snATAC_02_matrix/WT_PFC_W26_5kb.h5ad')

# %%

from ALLCools.clustering import remove_black_list_region, significant_pc_test, binarize_matrix, filter_regions, lsi


# %%
binarize_matrix(adata_atac_bin)

filter_regions(adata_atac_bin)


# %%
adata_atac_bin.var

# %%
split_rows = adata_atac_bin.var.index.str.split("-")
process_rows = []   
for i in split_rows:
    if len(i) == 3:
        process_rows.append([i[0],i[1],i[2]])
    else:
        print(i)
        process_rows.append(i[0]+i[1],i[2],i[3])


# %%
rows_p = np.array(process_rows)

# %%
adata_atac_bin.var['chrom'] = rows_p[:,0]
adata_atac_bin.var['start'] = rows_p[:,1].astype(int)
adata_atac_bin.var['end'] = rows_p[:,2].astype(int)


# %%
black_list_path

# %%
lsi(adata_atac_bin, algorithm='arpack', obsm='X_pca')



