import matplotlib.pyplot as plt
import anndata
import scanpy as sc
import scanpy.external as sce

from ALLCools.clustering import significant_pc_test
from ALLCools.plot import *
from harmonypy import run_harmony
import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(description="Calculate the methlyation level for a given region")
parser.add_argument("--input", type=str, required=False, help="Input path of the bw file folder",
                    default="atacsc-3region-clustering.h5ad")
parser.add_argument("--lamb", type=float, required=False, help="Lambda of harmony",default=0.3)
parser.add_argument("--nclust", type=int, required=False, help="Number of clusters",default=11)

args = parser.parse_args()
input_file = args.input

paramname = f"lamb{args.lamb}_nclust{args.nclust}"

adata_concat = anndata.read_h5ad(f'output/{input_file}')

default_n_threads = 1
os.environ['OPENBLAS_NUM_THREADS'] = f"{default_n_threads}"
os.environ['MKL_NUM_THREADS'] = f"{default_n_threads}"
os.environ['OMP_NUM_THREADS'] = f"{default_n_threads}"

# ho = run_harmony(adata_concat.obsm['X_pca'],
#                  meta_data=pd.DataFrame(adata_concat.obs['batch']),
#                  vars_use='batch',
#                  random_state=0,
#                  nclust=100,
#                  max_iter_harmony=20)

adata_sc = adata_concat[adata_concat.obs.sample_id.str.contains('sc')]
adata_atac = adata_concat[~adata_concat.obs.sample_id.str.contains('sc')]
sce.pp.mnn_correct(datas=[adata_sc, adata_atac], var_subset='highly_variable', batch_key='sample_id')
## lamb 11-12,lamb:0.3 is good
#adata_concat.obsm['X_pca'] = ho.Z_corr.T
# adata_concat.obs['umap_0'] = adata_concat.obsm['X_umap'][:, 0]
# adata_concat.obs['umap_1'] = adata_concat.obsm['X_umap'][:, 1]
# 

sc.pp.neighbors(adata_concat, n_neighbors=20,use_rep='X_pca_harmony')
sc.tl.leiden(adata_concat, resolution=1.5)
sc.tl.umap(adata_concat)

input_name = input_file.split('.h5ad')[0]

adata_concat.write_h5ad(f'output/{input_name}{paramname}.harmony.h5ad')

adata_concat.obs['umap_0'] = adata_concat.obsm['X_umap'][:, 0]
adata_concat.obs['umap_1'] = adata_concat.obsm['X_umap'][:, 1]
# 
fig, axes = plt.subplots(figsize=(8, 3), dpi=250, ncols=2)
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
                    hue='celltype.L1',
                    show_legend=True,
                    max_points=None,
                    s=1)

plt.savefig(f'output/{input_name}{paramname}.pdf')