import matplotlib.pyplot as plt
import anndata
import scanpy as sc

from ALLCools.clustering import significant_pc_test
from ALLCools.plot import *
from harmonypy import run_harmony
import pandas as pd
import os
import scanpy.external as sce
adata_concat = anndata.read_h5ad('output/merged-PFC-concat-clustering.h5ad')

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
sce.pp.harmony_integrate(adata_concat, 'batch')
#adata_concat.obsm['X_pca'] = ho.Z_corr.T
adata_concat.obs['umap_0'] = adata_concat.obsm['X_umap'][:, 0]
adata_concat.obs['umap_1'] = adata_concat.obsm['X_umap'][:, 1]
fig, axes = plt.subplots(figsize=(8, 3), dpi=250, ncols=2)

sc.pp.neighbors(adata_concat, n_neighbors=20,use_rep='X_pca_harmony')
sc.tl.leiden(adata_concat, resolution=1.5)
sc.tl.umap(adata_concat)

adata_concat.write_h5ad('output/merged-PFC-harmony-clustering.h5ad')
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

plt.savefig('output/merged-PFC-concat-clustering.pdf')