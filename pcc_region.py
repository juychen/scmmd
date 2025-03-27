# %%
import matplotlib.pyplot as plt
import anndata as adata
import scanpy as sc
import pandas as pd
import snapatac2 as snap
import numpy as np
from harmonypy import run_harmony
from pybedtools import BedTool
import subprocess
import argparse



parser = argparse.ArgumentParser(description="Calculate the methlyation level for a given region")
parser.add_argument("--region", type=str, required=False, help="Input path of the bw file folder",
                    default="HIP")
parser.add_argument("--celltype", type=str, required=False, help="Input path of the bw file folder",
                    default="Neuron")

args = parser.parse_args()

# %%
# regions = ['AMY','HIP','PFC']
# celltypes = ['OPC-Oligo', 'Immune','Astro-Epen','Vascular','Neuron']

# for celltype in celltypes:
#     for region in regions:
#         # read the network genertaed by CICRE
#         df_network = pd.read_csv(f'/data2st1/junyi/output/cicre/{region}_{celltype}_circe_network.csv', index_col=0)
#         distance = df_network.Peak1.str.split('_').str[1].astype(int) -df_network.Peak2.str.split('_').str[1].astype(int)
#         distance[distance.abs()<500_000].abs().hist(bins=100)
#         plt.axvline(x=np.median(distance[distance.abs()<500_000].abs()), color='r', linestyle='-' )  # y=position, color, linestyle
#        # plt.text(np.median(distance[distance.abs()<500_000].abs()) +10000, 25000, , color='k')  # Position text at (x=1, y=5.1)
#         plt.title(f'Median distance:{np.median(distance[distance.abs()<500_000].abs())}')
#         plt.savefig(f'/data2st1/junyi/output/cCRE/{region}_{celltype}_circe_network_distance.png')
#         plt.close()


# %%
# Read meta data and neighbor information
meta_sc = pd.read_csv('/data2st1/junyi/output/Brain.3regions.MetaData.csv.gz', index_col=0, compression='gzip')
meta_atac = pd.read_csv('/data2st1/junyi/output/ATAC.3regions.MetaData.csv.gz', index_col=0)
# Neighbor information
# Pseudo bulk to pseudo bulk level distances and mapping
df_nn = pd.read_csv('/data2st1/junyi/output/pseudo_sc_atac_nn.csv', index_col=None)

# %%
#/data2st1/junyi/output/cicre/HIP_OPC-Oligo_circe.h5ad
celltype = 'Neuron'
region = 'HIP'

celltype = args.celltype
region = args.region


# %%
# Read h5ad
adata_atac = sc.read_h5ad(f'/data2st1/junyi/output/cicre/{region}_{celltype}_circe.h5ad')

# %%
# Filter out the cells that are not in the meta data
adata_atac.obs['pseudo_group'] = meta_atac.loc[adata_atac.obs.index]['pseudo_group']

# %%
adata_atac.obs

# %%
# Aggregate the data
group_ata = snap.tl.aggregate_X(adata_atac,groupby='pseudo_group')

# %%
# Save the reagion information to bed file
adata_atac.var.iloc[:,1:4].to_csv(f'/data2st1/junyi/output/cCRE/atac_{region}_{celltype}.var.bed', sep='\t', header=False,index=False)

# %%
# Sort the bed file by chromosome and start position
command_str = f"sort -k1,1 -k2,2n /data2st1/junyi/output/cCRE/atac_{region}_{celltype}.var.bed > /data2st1/junyi/output/cCRE/atac_{region}_{celltype}_sort_var.bed"
subprocess.run(command_str, shell=True)


# %%
# Map the cCRE to the promoter
query = BedTool(f'/data2st1/junyi/output/cCRE/atac_{region}_{celltype}_sort_var.bed')
target2 = BedTool(f'/data2st1/junyi/generegion_vM33/promoter.sort.bed')
result = query.closest(target2, stream=True, D="ref")
result.saveas(f'/data2st1/junyi/output/cCRE/atac_{region}_{celltype}.promoter.bed')

# %%
# Read the promoter information
df_creG = pd.read_csv(f'/data2st1/junyi/output/cCRE/atac_{region}_{celltype}.promoter.bed', sep='\t', header=None)
df_creG.columns = ['chrom','start','end','chrom_p','start_p','end_p','score','strand','gene_name','gene_id','cCRE','dist']
df_creG['dist_up_p'] = df_creG['dist']
df_creG.loc[df_creG.strand=='-', 'dist_up_p'] = df_creG['dist']*-1

# %%
df_creG.loc[(df_creG.dist_up_p>0) & (df_creG.dist_up_p.abs()<500_000), 'dist_up_p'].hist(bins=100)
#plt.ylim(0,200_000)
#plt.xlim(0,500_000)


# %%


# %%
df_psuatac= df_nn.drop_duplicates('atac', keep='first')
df_psuatac['atac_pseudo'] = df_psuatac['atac_pseudo'].str.replace('-0','')
df_psuatac['sc_pseudo'] = df_psuatac['sc_pseudo'].str.replace('-1','')
group_ata.obs = df_psuatac.set_index('atac_pseudo').merge(group_ata.obs, left_index=True, right_index=True,how='right')

# %%
group_ata = group_ata[~group_ata.obs['sc_pseudo'].isna()]

# %%
df_cCREG = df_creG.drop_duplicates(['chrom','start','end'])

# %%
df_cCREG['region'] = df_cCREG['chrom'] + '_' + df_cCREG['start'].astype(str) + '_' + df_cCREG['end'].astype(str)
df_cCREG

# %%
group_ata.var = group_ata.var.merge(df_cCREG, left_index=True, right_on='region', how='left').set_index('region')

# %%
adata_sc_pse = sc.read_h5ad('/data2st1/junyi/output/PseudoCellAdata.3regions.h5ad')

# %%
adata_sc_pse.var.head()

# %%
group_ata.var.head()

# %%
group_ata

# %%
adata_scselected = adata_sc_pse[group_ata.obs['sc_pseudo']]

# %%
adata_scselected

# %%
from scipy.stats import pearsonr

rs = []
ps = []

count = 0
for idx, row in group_ata.var.iterrows():
    count += 1

    gid= row['gene_id'].split('.')[0]
    atac_X = group_ata[:,idx].X

    if gid in adata_scselected.var.index:
        sc_X = adata_sc_pse[group_ata.obs['sc_pseudo'],gid].X.todense()
    else:
        rs.append(np.nan)
        ps.append(1)
        continue
    
    if (np.count_nonzero(sc_X) == 0) or (np.count_nonzero(atac_X) == 0):
        rs.append(np.nan)
        ps.append(1)
        continue
    
    r,pva = pearsonr(atac_X, sc_X)
    rs.append(r[0])
    ps.append(pva[0])



# %%
pccs = np.nan_to_num(rs, nan=-1)
pvals =  np.nan_to_num(ps, nan=1)

# %%
idx_shullfed = group_ata.obs.sample(frac=1).index
group_ata_shuffled = group_ata[idx_shullfed]

# %%
rnull = []
pnull = []

count = 0
for idx, row in group_ata.var.iterrows():
    count += 1

    gid= row['gene_id'].split('.')[0]
    atac_X = group_ata_shuffled[:,idx].X
    if gid in adata_scselected.var.index:
        sc_X = adata_sc_pse[group_ata.obs['sc_pseudo'],gid].X.todense()
    else:
        # sc_X = np.zeros(atac_X.shape)
        rnull.append(np.nan)
        pnull.append(1)
        continue

    if (np.count_nonzero(sc_X) == 0) or (np.count_nonzero(atac_X) == 0):
        rnull.append(np.nan)
        pnull.append(1)
        continue

    r,pva = pearsonr(atac_X, sc_X)
    rnull.append(r[0])
    pnull.append(pva[0])

# %%
pccsnull = np.nan_to_num(rnull, nan=-1)
pvalsnull =  np.nan_to_num(pnull, nan=1)

# %%
df_pccs = pd.DataFrame({'pcc':pccs, 'pval':pvals, 'pcc_null':pccsnull, 'pval_null':pvalsnull}, index=group_ata.var.index)

# %%
df_pccs.to_csv(f'/data2st1/junyi/output/cCRE/{region}_{celltype}_circe_pccs.csv')

# %%
plt.hist(df_pccs['pcc_null'], bins=100, alpha=0.2, label='Null distribution')
plt.hist(df_pccs['pcc'], bins=100, alpha=0.5, label='PCC of gene and cCRE')
plt.legend()
plt.xlabel('PCC Null Values')
plt.ylabel('Frequency')
plt.title('Comparison of PCC Null Distributions')
plt.xlim(-0.5,0.5)
plt.axvline(np.percentile(df_pccs['pcc_null'],95), color='r', linestyle='-' )  # y=position, color, linestyle
plt.text(np.percentile(df_pccs['pcc_null'],95) +0.01, 2000, f'95% cutoff', color='r')  # Position text at (x=1, y=5.1)
plt.show()
plt.savefig(f'/data2st1/junyi/output/cCRE/{region}_{celltype}_circe_pccs.png')



