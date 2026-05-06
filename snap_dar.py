# %%
import matplotlib.pyplot as plt
import anndata
import scanpy as sc
import snapatac2 as snap
import numpy as np
import pandas as pd
import os
import scanpy.external as sce
#adata_concat = snap.read_dataset('output/motif/motif/motif/motif/motif/mouse_brain.h5ads')

# %%
file = "/data2st1/junyi/output/mouse_brain_dar.h5ad"

# %%
adata = anndata.read_h5ad(file)
region = 'ALL'
celltype = 'ALL'
base_name = f"{region}_{celltype}"
adata_AMY_neuron = adata

adata_AMY_neuron.obs['expriment'] = adata_AMY_neuron.obs['sample'].str[:2]


snap.tl.macs3(adata_AMY_neuron, groupby='celltype.L1.new')
peaks = snap.tl.merge_peaks(adata_AMY_neuron.uns['macs3'], snap.genome.GRCm39)

peak_mat = snap.pp.make_peak_matrix(adata_AMY_neuron, use_rep=peaks['Peaks'])

peak_mat.layers['raw'] = peak_mat.X.copy()
peak_mat.obsm['X_umap'] = adata_AMY_neuron.obsm['X_umap']
peak_mat.write(f"output/motif/{base_name}.h5ads")