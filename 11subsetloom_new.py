# %%
import anndata
import numpy as np
import os
import argparse
import loompy as lp

parser = argparse.ArgumentParser(
    description="Create cell-type-specific loom files from region-specific h5ad files"
)
parser.add_argument(
    "--input",
    nargs="+",
    required=True,
    metavar="REGION=H5AD",
    help=(
        "Region-specific h5ad inputs, for example "
        "PFC=/path/PFC.h5ad STR=/path/STR.h5ad. "
        "For one unlabelled path, use --region or the file name is used."
    ),
)
parser.add_argument("--output", type=str, required=False, help="Output path of the result",
                    default="/data2st2/junyi/output/stg1028/allsn/")
parser.add_argument(
    "--adata-output",
    type=str,
    default=None,
    help="Optional directory for the cell-type-specific h5ad files.",
)
parser.add_argument(
    "--region",
    type=str,
    default=None,
    help="Brain-region label for a single unlabelled --input path.",
)
parser.add_argument(
    "--celltype_column",
    type=str,
    required=False,
    help="Column for cell type",
    default="celltype.L2",
)


args = parser.parse_args()

def parse_input_specs(input_specs, fallback_region=None):
    parsed_specs = []
    for input_spec in input_specs:
        if "=" in input_spec:
            region_name, input_path = input_spec.split("=", 1)
            if not region_name or not input_path:
                raise ValueError(
                    f"Invalid --input value {input_spec!r}; use REGION=H5AD."
                )
        else:
            input_path = input_spec
            if len(input_specs) == 1 and fallback_region:
                region_name = fallback_region
            else:
                region_name = os.path.splitext(os.path.basename(input_path))[0]

        if region_name in {parsed_region for parsed_region, _ in parsed_specs}:
            raise ValueError(f"Duplicate brain-region label: {region_name}")
        parsed_specs.append((region_name, input_path))

    return parsed_specs


def build_celltype_region(obs_frame, region_name, celltype_column):
    if "region" in obs_frame.columns:
        region_values = obs_frame["region"].fillna(region_name).astype(str)
    else:
        region_values = region_name

    celltype_values = obs_frame[celltype_column].fillna("NA").astype(str)
    celltype_region = region_values + "_" + celltype_values
    celltype_region = celltype_region.str.replace('/', '-', regex=False)
    return celltype_region.str.replace(' ', '_', regex=False)


input_configs = parse_input_specs(args.input, args.region)
os.makedirs(args.output, exist_ok=True)
if args.adata_output:
    os.makedirs(args.adata_output, exist_ok=True)

datasets = []
all_celltypes = set()

try:
    for region_name, dataset_path in input_configs:
        print(f"Reading {region_name}: {dataset_path}")
        adata_backed = anndata.read_h5ad(dataset_path, backed="r")
        if args.celltype_column not in adata_backed.obs.columns:
            raise KeyError(
                f"{dataset_path} does not contain {args.celltype_column!r} in obs."
            )

        obs_frame = adata_backed.obs.copy()
        obs_frame["celltype.L2.region"] = build_celltype_region(
            obs_frame, region_name, args.celltype_column
        )
        all_celltypes.update(obs_frame["celltype.L2.region"].unique())
        datasets.append((region_name, adata_backed, obs_frame))

    for celltype in sorted(all_celltypes):
        try:
            print(f"Processing {celltype}")

            subset_list = []
            for _, adata_backed, obs_frame in datasets:
                selected_obs = obs_frame[obs_frame["celltype.L2.region"] == celltype]
                if selected_obs.empty:
                    continue
                adata_subset = adata_backed[selected_obs.index, :].to_memory()
                adata_subset.obs = selected_obs.copy()
                subset_list.append(adata_subset)

            if not subset_list:
                print(f"No cells found for {celltype}, skipping.")
                continue

            if len(subset_list) == 1:
                adata_subset = subset_list[0]
            else:
                adata_subset = anndata.concat(
                    subset_list, join="inner", merge="same", index_unique="-"
                )

            base_name = celltype.replace(" ", "_").replace("/", "-")
            loom_path = os.path.join(args.output, f"{base_name}.loom")

            if args.adata_output:
                adata_path = os.path.join(args.adata_output, f"{base_name}.h5ad")
                if os.path.exists(adata_path):
                    print(f"File {base_name}.h5ad already exists, skipping.")
                else:
                    adata_subset.write_h5ad(adata_path)

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
            print(f"Error processing {celltype}: {e}")
            continue
finally:
    for _, adata_backed, _ in datasets:
        adata_backed.file.close()
# %%

