import argparse
from pathlib import Path

import anndata as ad
import muon as mu
import numpy as np
import pandas as pd
import pychromvar as pc
import scanpy as sc
from pyjaspar import jaspardb


DEFAULT_GENOME = "/data2st1/junyi/ref/GRCm38.p6.genome.fa"
DEFAULT_JASPAR_RELEASE = "JASPAR2026"
DEFAULT_MIN_CELLS_PER_PEAK = 50
DEFAULT_MIN_GENES_PER_CELL = 2000
DEFAULT_MIN_TOTAL_COUNTS = 4000
DEFAULT_DOWNSAMPLE_CELLS = 10000


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run chromVAR deviations from an input h5ad file.")
    parser.add_argument("input_h5ad", help="Input h5ad file containing ATAC counts.")
    parser.add_argument("output_h5ad", help="Output h5ad file for chromVAR deviations.")
    parser.add_argument(
        "--genome-file",
        default=DEFAULT_GENOME,
        help=f"Reference genome FASTA path. Default: {DEFAULT_GENOME}",
    )
    parser.add_argument(
        "--jaspar-release",
        default=DEFAULT_JASPAR_RELEASE,
        help=f"JASPAR release to use. Default: {DEFAULT_JASPAR_RELEASE}",
    )
    parser.add_argument(
        "--min-cells-per-peak",
        type=int,
        default=DEFAULT_MIN_CELLS_PER_PEAK,
        help=f"Minimum cells per peak filter. Default: {DEFAULT_MIN_CELLS_PER_PEAK}",
    )
    parser.add_argument(
        "--min-genes-per-cell",
        type=int,
        default=DEFAULT_MIN_GENES_PER_CELL,
        help=f"Minimum genes-by-counts filter. Default: {DEFAULT_MIN_GENES_PER_CELL}",
    )
    parser.add_argument(
        "--min-total-counts",
        type=int,
        default=DEFAULT_MIN_TOTAL_COUNTS,
        help=f"Minimum total counts per cell. Default: {DEFAULT_MIN_TOTAL_COUNTS}",
    )
    parser.add_argument(
        "--no-downsample",
        action="store_true",
        help=f"Disable the default behavior of downsampling to at most {DEFAULT_DOWNSAMPLE_CELLS} cells after filtering.",
    )
    parser.add_argument(
        "--downsample-cells",
        type=int,
        default=DEFAULT_DOWNSAMPLE_CELLS,
        help=f"Maximum number of cells to keep when --downsample is enabled. Default: {DEFAULT_DOWNSAMPLE_CELLS}",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=0,
        help="Random seed used for downsampling. Default: 0",
    )
    return parser.parse_args()


def prepare_input_adata(adata: ad.AnnData, args: argparse.Namespace) -> ad.AnnData:
    if "count" not in adata.layers:
        raise KeyError("Input AnnData must contain a 'count' layer for chromVAR.")

    adata = adata.copy()
    adata.var["chr"] = adata.var_names.astype(str)
    adata.var_names = adata.var["chr"].str.replace(":", "-", regex=False)
    adata.X = adata.layers["count"].copy()

    sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
    mu.pp.filter_var(adata, "n_cells_by_counts", lambda values: values >= args.min_cells_per_peak)
    mu.pp.filter_obs(adata, "n_genes_by_counts", lambda values: values >= args.min_genes_per_cell)
    mu.pp.filter_obs(adata, "total_counts", lambda values: values >= args.min_total_counts)
    adata.X = adata.X.astype(np.float32)
    return adata


def maybe_downsample_adata(adata: ad.AnnData, args: argparse.Namespace) -> ad.AnnData:
    if args.no_downsample or adata.n_obs <= args.downsample_cells:
        return adata

    rng = np.random.default_rng(args.random_state)
    selected_obs = rng.choice(adata.obs_names.to_numpy(), size=args.downsample_cells, replace=False)
    selected_obs = pd.Index(selected_obs)
    return adata[selected_obs].copy()


def run_chromvar(adata: ad.AnnData, args: argparse.Namespace, output_path: Path) -> ad.AnnData:
    pc.add_peak_seq(adata, genome_file=args.genome_file)
    pc.add_gc_bias(adata)

    jdb_obj = jaspardb(release=args.jaspar_release)
    motifs = jdb_obj.fetch_motifs(collection="CORE", tax_group=["vertebrates"])

    pc.get_bg_peaks(adata)
    pc.match_motif(adata, motifs=motifs)

    df_motif_match = pd.DataFrame(
        adata.varm["motif_match"],
        index=adata.var_names,
        columns=adata.uns["motif_name"],
    )
    motif_match_path = output_path.parent / "motif_match.csv"
    df_motif_match.to_csv(motif_match_path)

    return pc.compute_deviations(adata)


def main() -> None:
    args = parse_args()
    input_path = Path(args.input_h5ad)
    output_path = Path(args.output_h5ad)

    print(f"Reading input AnnData: {input_path}")
    adata = sc.read_h5ad(input_path)
    print(f"Input shape before filtering: {adata.shape}")

    adata = prepare_input_adata(adata, args)
    print(f"Input shape after filtering: {adata.shape}")

    adata = maybe_downsample_adata(adata, args)
    print(f"Input shape for chromVAR: {adata.shape}")

    print("Running chromVAR...")
    dev = run_chromvar(adata, args, output_path)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Writing chromVAR output: {output_path}")
    dev.write_h5ad(output_path)
    print("Done.")


if __name__ == "__main__":
    main()
