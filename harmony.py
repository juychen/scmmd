import matplotlib.pyplot as plt
import anndata
import scanpy as sc

from ALLCools.clustering import significant_pc_test
from ALLCools.plot import *
from harmonypy import run_harmony
import pandas as pd
import os
import scanpy.external as sce

input_file = 'atacsc-3region-pseudo-clustering.h5ad'

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
sce.pp.harmony_integrate(adata_concat, 'batch',
                         **{'nclust':9,'lamb':1,'epsilon_harmony':0,'epsilon_cluster':0,'max_iter_harmony':10})
#adata_concat.obsm['X_pca'] = ho.Z_corr.T
# adata_concat.obs['umap_0'] = adata_concat.obsm['X_umap'][:, 0]
# adata_concat.obs['umap_1'] = adata_concat.obsm['X_umap'][:, 1]
# 

sc.pp.neighbors(adata_concat, n_neighbors=20,use_rep='X_pca_harmony')
sc.tl.leiden(adata_concat, resolution=1.5)
sc.tl.umap(adata_concat)

input_name = input_file.split('.h5ad')[0]

adata_concat.write_h5ad(f'output/{input_name}.harmony.h5ad')

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

plt.savefig(f'output/{input_name}.pdf')