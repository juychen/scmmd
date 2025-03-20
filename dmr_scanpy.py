# %%
# import dask.dataframe as dd


# # %%
# ddf = dd.read_csv('/data2st1/junyi/snmc_bedgraph/merged/all_unified.bedGraph',sep='\t',header=0, assume_missing=True)

# # %%
# ddf

# # %%
# matrix = ddf.iloc[:,3:]

# # %%
# matrix = matrix.compute()

# # %%
# import scipy.sparse as sp
# sparse_matrix = sp.csr_matrix(matrix.values)


# # %%
# import scanpy as sc
# adata = sc.AnnData(sparse_matrix)

# # %%
# adata = adata.T

# # %%
# region_df = ddf.iloc[:,0:3]

# # %%
# region_df = region_df.compute()

# # %%
# region_df['start'] = region_df['start'].astype(int)
# region_df['end'] = region_df['end'].astype(int)

# # %%
# region_df

# # %%
# adata.obs_names = ddf.columns[3:]

# # %%
# adata.var = region_df

# # %%
# adata.var.index = adata.var['chrom'] + ':' + adata.var['start'].astype(str) + '-' + adata.var['end'].astype(str)

# # %%
# adata.var.index = adata.var.index.astype(str)
# adata.obs.index = adata.obs.index.astype(str)

# # %%
# adata.write_h5ad('/data2st1/junyi/snmc_bedgraph/merged/meth_frac_all.h5ad')
import scanpy as sc

adata = sc.read_h5ad('/data2st1/junyi/snmc_bedgraph/merged/meth_frac_all.h5ad')

# %%
mapping = {
    'Oligo_NN': 'OPC-Oligo',
    'Microglia_NN': 'Immune',
    'Endo_NN': 'Vascular',
    'Astro-TE_NN': 'Astro-Epen',
    'Astro-OLF_NN': 'Astro-Epen',
    'Astro-NT_NN': 'Astro-Epen',
    'Astroependymal_NN': 'Astro-Epen',
    'Hypendymal_NN': 'Astro-Epen',
    'Ependymal_NN': 'Astro-Epen',
    'Tanycyte_NN': 'Astro-Epen',
    'Bergmann_NN': 'Astro-Epen',
    'CHOR_NN': 'Astro-Epen',
    'Astro-CB_NN': 'Astro-Epen',
    'OPC_NN': 'OPC-Oligo',
    'OEC_NN': 'Astro-Epen',
    'BAM_NN': 'Immune',
    'DC_NN': 'Immune',
    'Lymphoid_NN': 'Immune',
    'VLMC_NN': 'Vascular',
    'ABC_NN': 'Immune',
    'SMC_NN': 'Vascular',
    'Peri_NN': 'Vascular'
}

# %%
adata.obs['Subtype'] = adata.obs.index.str.replace('.CGN','')

# %%
adata.obs['Celltype.l1'] = adata.obs['Subtype'].map(mapping).fillna('Neuron')


# %%
adata.obs['Celltype.l2'] = adata.obs['Subtype'].str.split('_').str[0]+"_"+adata.obs['Subtype'].str.split('_').str[-1]

# %%
#sc.pp.highly_variable_genes(adata, n_top_genes=5000)
sc.tl.rank_genes_groups(adata, 'Celltype.l1', method='wilcoxon',pts=True)

celltypes = ['OPC-Oligo', 'Immune','Astro-Epen','Vascular','Neuron']

for celltype in celltypes:
    df = sc.get.rank_genes_groups_df(adata, group=celltype, key='rank_genes_groups',pval_cutoff=0.05)
    df.to_csv(f"/data2st1/junyi/snmc_bedgraph/merged/DMR_{celltype}_wilcoxon.csv")
    #sc.pl.umap(adata, color=df.sort_values('logfoldchanges',ascending=False).head(10).names, size=50,save=f"ALL_{celltype}_umap.png")

adata.write_h5ad('/data2st1/junyi/snmc_bedgraph/merged/meth_frac_all.h5ad')
