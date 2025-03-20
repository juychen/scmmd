# %%
import numpy as np
import circe as ci
import scanpy as sc
import scipy as sp
import warnings
import os
warnings.filterwarnings('ignore')

# %%
regions = ['AMY','HIP','PFC']
celltypes = ['OPC-Oligo', 'Immune','Astro-Epen','Vascular','Neuron']

for region in regions:
    for celltype in celltypes:
        # %%
        # region = 'AMY'
        # celltype = 'Vascular'
        base_name = f"{region}_{celltype}"

        if os.path.exists(f"/data2st1/junyi/output/cicre/{base_name}_circe.h5ad"):
            continue



        adata = sc.read_h5ad(f"/data2st1/junyi/output/motif/{base_name}.h5ads")

        # %%
        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.filter_cells(adata, min_genes=200)


        # %%
        adata.var.index=adata.var.index.str.replace(":", "_").str.replace("-", "_")

        # %%
        # Column to use for stratification
        if adata.shape[0] > 10000:
            stratify_column = 'sample'
            stratify_column = 'sample'

            # Number of cells to sample from each group
            n_samples_per_group = 2000

            # Perform stratified sampling
            sampled_indices = (
                adata.obs
                .groupby(stratify_column, group_keys=False)
                .apply(lambda x: x.sample(min(n_samples_per_group, len(x))))
                .index
            )
            adata  = adata[sampled_indices, :]


        adata = ci.add_region_infos(adata)

        ci.compute_atac_network(
            adata, #metacells,
            organism="mouse",
        )

        final_score = ci.sliding_graphical_lasso(
            adata,
            n_samples=50,
            n_samples_maxtry=100,
            max_alpha_iteration=500,
            verbose=True
        )
        adata.varp['atac_network'] = final_score


        # %%
        circe_network = ci.extract_atac_links(adata) #metacells)


        # %%
        circe_network.to_csv(f"/data2st1/junyi/output/cicre/{base_name}_circe_network.csv")

        # %%
        ccans = ci.find_ccans(circe_network, seed=0)


        # %%
        ccans.to_csv(f"/data2st1/junyi/output/cicre/{base_name}_ccans.csv")

        # %%
        adata = ci.add_ccans(adata)
        
        try:
            adata.var['CCAN'] = adata.var['CCAN'].astype(str)
            adata.write(f"/data2st1/junyi/output/cicre/{base_name}_circe.h5ad")
        except:
            print(f"Error in wirting h5ad {base_name}")
            continue



