# %%
import matplotlib.pyplot as plt
import snapatac2 as snap
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import silhouette_score

import warnings
warnings.filterwarnings("ignore")


# %%
adata_concat = snap.read_dataset('/data2st1/junyi/output/atac0416/mouse_brain_subset.h5ads')

# %%
df_meta_all = pd.read_csv('/data2st1/junyi/output/atac0416/ATACSC_3REGION_ALL_L2annoated.csv',index_col=0)

# %%
black_list = ['Immune','OPC-Oligo','Doublet','PFC Doublet','PFC Not sure','Not sure','Astro-Epen']

# %%
# df_meta_all["celltype.L2"].unique()
# l2type = df_meta_all["celltype.L2"].unique()[0]

# %%
best_ress = []
n_clusters = []
n_cells = []
ctypes = []

for l2type in df_meta_all["celltype.L2"].unique():
    if l2type in black_list:
        continue

    # Subset the data for the current L2 type
    st1 = df_meta_all[df_meta_all["celltype.L2"]==l2type]
    l2name = l2type.replace(' ','_')
    datasubset,newidx =adata_concat.subset(obs_indices=st1.index,out=f'/data2st1/junyi/output/atac0416/iterative/l3/{l2name}.h5ads')

    # Preprocess the subset 
    snap.tl.spectral(datasubset,n_comps=30)
    snap.pp.knn(datasubset, use_rep="X_spectral",n_neighbors=50)
    snap.tl.umap(datasubset)

    # Perform clustering with different resolutions and calculate silhouette scores
    max_silhoutte =0
    X=datasubset.obsm['X_spectral']
    df_plot = pd.DataFrame()

    # Iterate over different resolution values
    for res in np.arange(0.1, 2.0, 0.1):
        snap.tl.leiden(datasubset, resolution=res)
        Y=datasubset.obs['leiden']

        # Calculate silhouette score
        try:
            score = silhouette_score(X,Y)
        except ValueError as e:
            score = 0
        if score > max_silhoutte:
            max_silhoutte = score
            best_res = res
            print(f"Best silhouette score: {max_silhoutte} at resolution {best_res}")

        if len(df_plot) == 0:
            df_plot = pd.DataFrame({"UMAP1":datasubset.obsm['X_umap'][:,0],"UMAP2":datasubset.obsm['X_umap'][:,1],
                                    "resolution":np.round(res,1),"silhouette_score":score,
                                    "leiden":Y})

        else:
            df_tmp = pd.DataFrame({"UMAP1":datasubset.obsm['X_umap'][:,0],"UMAP2":datasubset.obsm['X_umap'][:,1],
                                    "resolution":np.round(res,1),"silhouette_score":score,
                                    "leiden":Y})
            df_plot = pd.concat([df_plot,df_tmp],ignore_index=True)
        
        # Stop if the median cluster size is less than 400
        if Y.value_counts().to_pandas()["count"].median()<400:
            print(f"Resolution {res} has a median cluster with less than 400 cells")
            break

    # Plot the results
    snap.tl.leiden(datasubset, resolution=best_res)
    sns.relplot(
        data=df_plot, x="UMAP1", y="UMAP2",
        col="resolution", hue="leiden",
        col_wrap=5,linewidth=0,
        height=3, aspect=1.3
    )
    plt.savefig(f'/data2st1/junyi/output/atac0416/iterative/l3/{l2name}_umap.png', dpi=300, bbox_inches='tight')
    df_meta_subtypes = pd.DataFrame({"celltype.L3":datasubset.obs['leiden'].to_numpy()},index=datasubset.obs_names)
    df_meta_subtypes["celltype.L3"] = l2type+"-"+df_meta_subtypes["celltype.L3"].astype(str)
    df_meta_subtypes.to_csv(f'/data2st1/junyi/output/atac0416/iterative/l3/{l2name}_l3.csv')

    best_ress.append(best_res)
    n_clusters.append(len(datasubset.obs['leiden'].unique()))
    ctypes.append(l2type)
    n_cells.append(len(datasubset.obs_names))
    datasubset.close()

# %%
df_summary = pd.DataFrame({"celltype.L2":ctypes,"best_resolution":best_ress,"n_clusters":n_clusters,"n_cells":n_cells})  

# %%
df_summary.to_csv('/data2st1/junyi/output/atac0416/iterative/l3/l3_summary.csv')


