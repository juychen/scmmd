# %%
import anndata
import numpy as np
import os
import argparse
import loompy as lp

parser = argparse.ArgumentParser(
    description=(
        "Create per-celltype loom files from multiple 3v3 downsampled h5ad files, "
        "keeping all conditions and de-duplicating shared control-group cells."
    )
)
parser.add_argument(
    "--input",
    type=str,
    nargs="+",
    default=[
        "/data7/mark/STG/dataset/snRNA/merge_SCH_new/CUSUS_3v3_500_1000gene_new_final/adata_SUS_3v3_downsampled_ratio.h5ad",
        "/data7/mark/STG/dataset/snRNA/merge_SCH_new/CURES_3v3_500_1000gene_new/adata_RES_3v3_downsampled_ratio.h5ad",
        "/data7/mark/STG/dataset/snRNA/merge_SCH_new/CSSUS_3v3_500_1000gene_new/adata_CSDS_3v3_downsampled_ratio.h5ad",
        "/data7/mark/STG/dataset/snRNA/merge_SCH_new/CSRES_3v3_500_1000gene_new/adata_CSRES_3v3_downsampled_ratio.h5ad",
    ],
    help="List of region-specific h5ad inputs (REGION=/path/file.h5ad or plain path).",
)
parser.add_argument(
    "--output", type=str,
    default="/data2st2/junyi/output/stg1028/allsn_3v3/",
    help="Output directory for loom files.",
)
parser.add_argument(
    "--adata-output", type=str, default=None,
    help="Optional directory for merged per-celltype h5ad files.",
)
parser.add_argument(
    "--celltype_column", type=str, default="celltype.L2",
    help="obs column holding the cell type annotation.",
)
parser.add_argument(
    "--dedup_key", type=str, default=None,
    help="obs column used to identify duplicate control cells. "
         "Default: use obs_names (cell barcodes). Set to e.g. 'cell_id' if barcodes differ.",
)

args = parser.parse_args()


def pd_index_to_str(idx):
    return idx.astype(str)


def parse_input_spec(input_spec):
    if "=" in input_spec:
        region_name, input_path = input_spec.split("=", 1)
        if not region_name or not input_path:
            raise ValueError(f"Invalid --input value {input_spec!r}; use REGION=H5AD.")
    else:
        input_path = input_spec
        region_name = os.path.splitext(os.path.basename(input_path))[0]
    return region_name, input_path


def build_target_names(obs_frame, region_name, celltype_column):
    if "region" in obs_frame.columns:
        region_values = obs_frame["region"].fillna(region_name).astype(str)
    else:
        region_values = region_name
    target = region_values + "_" + obs_frame[celltype_column].fillna("NA").astype(str)
    return target.str.replace("/", "-", regex=False).str.replace(" ", "_", regex=False)


os.makedirs(args.output, exist_ok=True)
if args.adata_output:
    os.makedirs(args.adata_output, exist_ok=True)

datasets = []          # (name, adata_backed, obs_frame)
all_targets = set()
seen_keys = None       # running set of dedup keys to drop repeated controls

try:
    # ---- Pass 1: read obs of every dataset, build target labels, mark duplicates ----
    for spec in args.input:
        name, path = parse_input_spec(spec)
        print(f"Reading {name}: {path}")
        adata_backed = anndata.read_h5ad(path, backed="r")
        if args.celltype_column not in adata_backed.obs.columns:
            raise KeyError(
                f"{path} does not contain {args.celltype_column!r} in obs. "
                f"Available: {list(adata_backed.obs.columns)}"
            )

        obs_frame = adata_backed.obs.copy()
        keys = (
            obs_frame[args.dedup_key].astype(str)
            if args.dedup_key is not None
            else pd_index_to_str(obs_frame.index)
        )
        dup_mask = keys.isin(seen_keys) if seen_keys is not None else np.zeros(len(obs_frame), dtype=bool)
        n_dup = int(dup_mask.sum())
        if n_dup:
            print(f"  {name}: dropping {n_dup} duplicated (control) cells.")
        obs_frame = obs_frame[~dup_mask].copy()
        # remember which cells were kept so later datasets drop them again
        kept_keys = keys[~dup_mask]
        seen_keys = set(kept_keys) if seen_keys is None else seen_keys | set(kept_keys)

        obs_frame["target.celltype"] = build_target_names(obs_frame, name, args.celltype_column)
        all_targets.update(obs_frame["target.celltype"].unique())
        datasets.append((name, adata_backed, obs_frame))

    # ---- Pass 2: subset each celltype across all datasets and write loom ----
    for target in sorted(all_targets):
        try:
            print(f"Processing {target}")
            subset_list = []
            for name, adata_backed, obs_frame in datasets:
                selected_obs = obs_frame[obs_frame["target.celltype"] == target]
                if selected_obs.empty:
                    continue
                sub = adata_backed[selected_obs.index, :].to_memory()
                sub.obs = selected_obs.copy()
                subset_list.append(sub)

            if not subset_list:
                print(f"No cells found for {target}, skipping.")
                continue

            adata_subset = (
                subset_list[0] if len(subset_list) == 1
                else anndata.concat(subset_list, join="inner", merge="same")
            )
            print(f"  {target}: {adata_subset.n_obs} cells from {len(subset_list)} dataset(s)")

            base_name = target.replace(" ", "_").replace("/", "-")
            loom_path = os.path.join(args.output, f"{base_name}.loom")

            if args.adata_output:
                adata_path = os.path.join(args.adata_output, f"{base_name}.h5ad")
                if not os.path.exists(adata_path):
                    adata_subset.write_h5ad(adata_path)
                else:
                    print(f"File {base_name}.h5ad already exists, skipping.")

            if os.path.exists(loom_path):
                print(f"File {base_name}.loom already exists, skipping.")
                continue

            matrix = adata_subset.X
            row_attrs = {"Gene": np.asarray(adata_subset.var_names)}
            col_attrs = {
                "CellID": np.asarray(adata_subset.obs_names),
                "nGene": np.asarray((matrix > 0).sum(axis=1)).flatten(),
                "nUMI": np.asarray(matrix.sum(axis=1)).flatten(),
            }
            lp.create(loom_path, matrix.transpose(), row_attrs, col_attrs)

        except Exception as e:
            print(f"Error processing {target}: {e}")
            continue

finally:
    for _, adata_backed, _ in datasets:
        try:
            adata_backed.file.close()
        except Exception:
            pass

# %%
