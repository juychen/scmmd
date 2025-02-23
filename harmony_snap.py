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
snap.pp.harmony(adata_concat, batch="sample", max_iter_harmony=20)

adata_concat.close()