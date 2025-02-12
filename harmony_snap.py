import matplotlib.pyplot as plt
import anndata
import scanpy as sc
import snapatac2 as snap
import numpy as np
import pandas as pd
import os
import scanpy.external as sce
adata_concat = snap.read_dataset('output/mouse_brain.h5ads')

default_n_threads = 10
os.environ['OPENBLAS_NUM_THREADS'] = f"{default_n_threads}"
os.environ['MKL_NUM_THREADS'] = f"{default_n_threads}"
os.environ['OMP_NUM_THREADS'] = f"{default_n_threads}"

# ho = run_harmony(adata_concat.obsm['X_pca'],
#                  meta_data=pd.DataFrame(adata_concat.obs['batch']),
#                  vars_use='batch',
#                  random_state=0,
#                  nclust=100,
#                  max_iter_harmony=20)
snap.pp.mnc_correct(adata_concat, batch="sample")
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