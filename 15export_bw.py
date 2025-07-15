# %%
import snapatac2 as snap
import scanpy as sc


# %%

# adata = sc.read_h5ad("/data2st1/junyi/output/atac0627/3REGIONS_peak.h5ads")
# adata.obs['celltype.L2.Condition'] = adata.obs['celltype.L2'].astype(str) + "_" + adata.obs['Condition'].astype(str)
# adata.write_h5ad("/data2st1/junyi/output/atac0627/3REGIONS_peak.h5ads")
# adata.obs.to_csv("/data2st1/junyi/output/atac0627/3REGIONS_peak_final.obs.csv")



# %%


# %%
data_frag = snap.read_dataset("/data2st1/junyi/output/atac0627/doublet_filtered.h5ads/_dataset.h5ads")

snap.ex.export_coverage(data_frag, groupby='celltype.L2.Condition',out_dir='/data2st1/junyi/output/atac0627/tracks', suffix='.bw')

# %%
data_frag.close()