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
    required=True,
    metavar="REGION=H5AD",
    help=(
        "One region-specific h5ad input, for example "
        "PFC=/path/PFC.h5ad. For an unlabelled path, use --region."
    ),
)
parser.add_argument("--output", type=str, required=False, help="Output path of the result",
                    default="/data2st2/junyi/output/stg1028/scvisn/")
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
parser.add_argument(
    "--layer",
    type=str,
    default=None,
    help="Expression layer to use; omitted means use adata.X.",
)


args = parser.parse_args()

def parse_input_spec(input_spec, fallback_region=None):
    if "=" in input_spec:
        region_name, input_path = input_spec.split("=", 1)
        if not region_name or not input_path:
            raise ValueError(
                f"Invalid --input value {input_spec!r}; use REGION=H5AD."
            )
    else:
        input_path = input_spec
        region_name = fallback_region or os.path.splitext(
            os.path.basename(input_path)
        )[0]

    return region_name, input_path


def build_celltype_region(obs_frame, region_name, celltype_column):
    if "region" in obs_frame.columns:
        region_values = obs_frame["region"].fillna(region_name).astype(str)
    else:
        region_values = region_name

    celltype_values = obs_frame[celltype_column].fillna("NA").astype(str)
    celltype_region = region_values + "_" + celltype_values
    celltype_region = celltype_region.str.replace('/', '-', regex=False)
    return celltype_region.str.replace(' ', '_', regex=False)


region_name, dataset_path = parse_input_spec(args.input, args.region)
os.makedirs(args.output, exist_ok=True)
if args.adata_output:
    os.makedirs(args.adata_output, exist_ok=True)

try:
    print(f"Reading {region_name}: {dataset_path}")
    adata_backed = anndata.read_h5ad(dataset_path, backed="r")
    if args.celltype_column not in adata_backed.obs.columns:
        raise KeyError(
            f"{dataset_path} does not contain {args.celltype_column!r} in obs."
        )
    if args.layer is not None and args.layer not in adata_backed.layers:
        raise KeyError(
            f"{dataset_path} does not contain layer {args.layer!r}. "
            f"Available layers: {list(adata_backed.layers.keys())}"
        )

    obs_frame = adata_backed.obs.copy()
    obs_frame["celltype.L0"] = np.where(
        obs_frame["Neurotransmitter"].eq("NN"),
        "NN",
        "N",
    )
    obs_frame["celltype.L0"] = obs_frame["celltype.L0"].astype("category")
    obs_frame["target.celltype"] = build_celltype_region(
        obs_frame, region_name, args.celltype_column
    )

    for celltype in sorted(obs_frame["target.celltype"].unique()):
        try:
            print(f"Processing {celltype}")

            selected_obs = obs_frame[obs_frame["target.celltype"] == celltype]
            if selected_obs.empty:
                print(f"No cells found for {celltype}, skipping.")
                continue

            adata_subset = adata_backed[selected_obs.index, :].to_memory()
            adata_subset.obs = selected_obs.copy()

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

            if args.layer is None:
                matrix = adata_subset.X
            else:
                matrix = adata_subset.layers[args.layer]
            
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
            if "adata_backed" in locals():
                adata_backed.file.close()
# %%

