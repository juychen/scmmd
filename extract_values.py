# %%
import pandas as pd
import tarfile
import io
from ALLCools.mcds import MCDS
import pyBigWig
import pandas as pd# %%
from wmb import cemba, mm10
import glob as glob
import numpy as np
import os
import matplotlib.pyplot as plt

# %%
black_list_path = mm10.ENCODE_BLACKLIST_PATH


# %%
black_list_path

# %%
mch_pattern = 'CHN'
mcg_pattern = 'CGN'
n_top_feature = 20000

# %%
gene_meta_path = '/home/junyichen/code/whole_mouse_brain/wmb/files/modified_gencode.vM23.primary_assembly.annotation.gene.flat.tsv.gz'

# %%
chrom_to_remove = ['chrM', 'chrX', 'chrY']


mcds_expression = MCDS.open(
    '/data2st1/junyi/methlyatlas/mCseq/data.nemoarchive.org/biccn/grant/u19_cemba/ecker/epigenome/cellgroup/mCseq3/mouse/processed/counts/CEMBA.snmC.L4Region.AIBS_TENX.log1pCPM.zarr/',
     obs_dim='CellGroup'
)


# %%
var_dim = "geneslop2k-vm23"
mcds = MCDS.open('/data2st1/junyi/methlyatlas/mCseq/data.nemoarchive.org/biccn/grant/u19_cemba/ecker/epigenome/cellgroup/mCseq3/mouse/processed/counts/CEMBA.snmC.L4RegionAgg.zarr',
                   obs_dim='L4Region',var_dim=var_dim)


# %%
abs(mcds['geneslop2k-vm23_start']-mcds['geneslop2k-vm23_end']).min()

# %%
use_gene_meta = pd.read_csv("/home/junyichen/code/scmmd/output/GeneMetadata.csv.gz",index_col=0)

# %%# %%
df_region = pd.read_csv('/data2st1/junyi/methlyatlas/atac/Subclass.names.map.csv.gz',index_col=0)

# %%
df_region

# %%
df_meta_integrateion = pd.read_csv("/data2st1/junyi/methlyatlas/mCseq/CEMBA.mC.Metadata/CEMBA.mC.Metadata.csv")

# %%
df_meta_integrateion.drop_duplicates(subset='CellGroup',inplace=True)

# %%
df_meta_integrateion.set_index('CellGroup',inplace=True)

# %%
df_meta_integrateion.loc[mcds["L4Region"]]

# %%
df_meta_integrateion


bw_files_mseq = glob.glob("/data2st1/junyi/methlyatlas/mCseq/data.nemoarchive.org/biccn/grant/u19_cemba/ecker/epigenome/cellgroup/mCseq3/mouse/processed/other/*CGN*bw")

cell_groups = []
methtypes =  []

for bw_file in bw_files_mseq:
    cell_group = os.path.basename(bw_file).split('.')[0]
    cell_groups.append(cell_group)
    methtype = os.path.basename(bw_file).split('.')[1]
    methtypes.append(methtype)
    bw = pyBigWig.open(bw_file)

    # %%
    bw.isBigWig()

    # %%
    use_gene_meta

    # %%
    stats =[]
    CGN_values = []
    for gene in use_gene_meta.index:
        chrom = use_gene_meta.loc[gene, 'chrom']
        start = use_gene_meta.loc[gene, 'start']
        end = use_gene_meta.loc[gene, 'end']
        
        values = bw.values(chrom, start-2000, start)
        arr_values = np.array(values)
        arr_values = np.nan_to_num(arr_values)
        stats.append(arr_values.mean())
        CGN_values.append(values)


    if not os.path.exists(f'output/{cell_group}'):
        os.makedirs(f'output/{cell_group}')
    plt.hist(stats, bins=100, color='skyblue', edgecolor='black')
    # Add labels and title
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.title('Histogram of Random Data')
    plt.savefig(f'output/{cell_group}/{methtype}.png')
    plt.close()

    value_array = np.array(CGN_values)
    np.save(f'output/{cell_group}/{methtype}.npy',value_array)
    stats = np.array(stats)
    np.save(f'output/{cell_group}/{methtype}_stats.npy',stats)



bw_files_atac = glob.glob("/data2st1/junyi/methlyatlas/atac/catlas.org/renlab_downloads/wholemousebrain/bigwig_macs2/*bw")

cell_groups = []
methtypes =  []

for bw_file in bw_files_atac:
    cell_group = os.path.basename(bw_file).split('.')[0]
    cell_groups.append(cell_group)
    bw = pyBigWig.open(bw_file)

    # %%
    bw.isBigWig()

    stats =[]
    atacpeak_values = []
    for gene in use_gene_meta.index:
        chrom = use_gene_meta.loc[gene, 'chrom']
        start = use_gene_meta.loc[gene, 'start']
        end = use_gene_meta.loc[gene, 'end']
        
        values = bw.values(chrom, start-2000, start)
        arr_values = np.array(values)
        arr_values = np.nan_to_num(arr_values)
        stats.append(arr_values.sum())
        atacpeak_values.append(values)


    if not os.path.exists(f'output/{cell_group}'):
        os.makedirs(f'output/{cell_group}')
    plt.hist(stats, bins=100, color='skyblue', edgecolor='black')
    # Add labels and title
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.ylim(0,1000)
    plt.title('Histogram of Random Data')
    plt.savefig(f'output/{cell_group}/atac.png')
    plt.close()

    value_array = np.array(CGN_values)
    np.save(f'output/{cell_group}/atac.npy',value_array)
    stats = np.array(stats)
    np.save(f'output/{cell_group}/atac_stats.npy',stats)