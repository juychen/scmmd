import warnings
warnings.filterwarnings("ignore")

import anndata as ad
import pandas as pd
import scanpy as sc
import scvi as sv
import numpy as np
import os
from matplotlib import pyplot as plt

input_file = 'atacsc-3region-pseudo-clustering.h5ad'

adata_concat = sc.read_h5ad(f'output/{input_file}')

default_n_threads = 1
os.environ['OPENBLAS_NUM_THREADS'] = f"{default_n_threads}"
os.environ['MKL_NUM_THREADS'] = f"{default_n_threads}"
os.environ['OMP_NUM_THREADS'] = f"{default_n_threads}"

# ho = run_harmony(adata_concat.obsm['X_pca'],
#                  meta_data=pd.DataFrame(adata_concat.obs['batch']),
#                  vars_use='batch',
#                  random_state=0,
#                  nclust=100,
# #                  max_iter_harmony=20)
# sce.pp.harmony_integrate(adata_concat, 'batch',
#                          **{'nclust':8,'lamb':0.5,'epsilon_harmony':0,'epsilon_cluster':0,'max_iter_harmony':10})

sv.model.SCVI.setup_anndata(adata_concat, batch_key="batch")
vae = sv.model.SCVI(
    adata_concat,
    n_layers=2,
    n_latent=30,
    gene_likelihood="nb",
    dispersion="gene-batch",
)

vae.train(max_epochs=1000, early_stopping=True)
ax = vae.history['elbo_train'][1:].plot()
vae.history['elbo_validation'].plot(ax=ax)
plt.savefig(f'output/{input_file}.scvi_training.pdf')

#adata_concat.obsm['X_pca'] = ho.Z_corr.T
# adata_concat.obs['umap_0'] = adata_concat.obsm['X_umap'][:, 0]
# adata_concat.obs['umap_1'] = adata_concat.obsm['X_umap'][:, 1]
# 
adata_concat.obs["celltype_scanvi"] = 'Unknown'
ref_idx = adata_concat.obs['batch'] == "sc"
adata_concat.obs["celltype_scanvi"][ref_idx] = adata_concat.obs['celltype.L1'][ref_idx]

lvae = sv.model.SCANVI.from_scvi_model(
    vae,
    adata=adata_concat,
    labels_key="celltype_scanvi",
    unlabeled_category="Unknown",
)
lvae.train(max_epochs=1000, n_samples_per_label=100)

adata_concat.obs["C_scANVI"] = lvae.predict(adata_concat)
adata_concat.obsm["X_scANVI"] = lvae.get_latent_representation(adata_concat)
sc.pp.neighbors(adata_concat, n_neighbors=20,use_rep='X_scANVI')
sc.tl.leiden(adata_concat, resolution=1.5)
sc.tl.umap(adata_concat)

input_name = input_file.split('.h5ad')[0]

adata_concat.write_h5ad(f'output/{input_name}.scanvi.h5ad')

adata_concat.obs['umap_0'] = adata_concat.obsm['X_umap'][:, 0]
adata_concat.obs['umap_1'] = adata_concat.obsm['X_umap'][:, 1]


sc.pl.umap(adata_concat, color=['batch','celltype.L1','C_scANVI'], save=f'{input_name}scanvi.pdf')