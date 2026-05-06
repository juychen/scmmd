#!/usr/bin/env python3
"""
Compute peak-gene Pearson correlations across pseudobulk samples.

Usage:
  python peak_gene_corr_pseudobulk.py \
    --adata-bulk /path/to/adata_bulk.h5ad \
    --adata-atac /path/to/adata_atac.h5ad \
    --promoter-peak /path/to/promoter_peak_500kb.csv \
    --out-prefix /path/to/output_prefix \
    --chunks 20

The script loads full gene list (no DEG filtering), builds candidate peak-gene pairs
from the promoter-peak CSV, computes Pearson correlations across shared obs (samples),
applies BH FDR correction, and writes both the full correlation table and the
significant subset (FDR < 0.05) as CSV files.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from scipy import sparse
from statsmodels.stats.multitest import multipletests
import scanpy as sc


def parse_args():
    p = argparse.ArgumentParser(description="Peak-Gene correlation on pseudobulk samples")
    p.add_argument('--adata-bulk', required=True, help='h5ad file for pseudobulk RNA')
    p.add_argument('--adata-atac', required=True, help='h5ad file for pseudobulk ATAC')
    p.add_argument('--promoter-peak', required=True, help='CSV from bedtools promoter-peak (500kb)')
    p.add_argument('--out-prefix', required=True, help='Output prefix for result files')
    p.add_argument('--chunks', type=int, default=20, help='Number of chunks/threads (default 20)')
    p.add_argument('--min-samples', type=int, default=3, help='Minimum samples required to compute correlation')
    return p.parse_args()


def to_dense(x):
    if sparse.issparse(x):
        return np.asarray(x.toarray())
    return np.asarray(x)


def process_chunk(chunk_df, peak_value_map, gene_value_map, min_samples=3):
    results = []
    for row in chunk_df.itertuples(index=False):
        peak = row.peak_id
        gene = row.gene_name
        pv = peak_value_map.get(peak)
        gv = gene_value_map.get(gene)
        if pv is None or gv is None:
            continue
        pv = np.asarray(pv).ravel()
        gv = np.asarray(gv).ravel()
        mask = np.isfinite(pv) & np.isfinite(gv)
        pv2 = pv[mask]
        gv2 = gv[mask]
        if pv2.size < min_samples:
            continue
        if np.all(pv2 == pv2[0]) or np.all(gv2 == gv2[0]):
            continue
        try:
            r, pval = pearsonr(pv2, gv2)
        except Exception:
            continue
        results.append({'peak_id': peak, 'gene_name': gene, 'correlation': float(r), 'pval': float(pval), 'n_samples': int(pv2.size)})
    return results


def main():
    args = parse_args()
    out_prefix = Path(args.out_prefix)

    print('Loading data...')
    adata_rna = sc.read_h5ad(args.adata_bulk)
    adata_atac = sc.read_h5ad(args.adata_atac)

    print('Loading promoter-peak table...')
    df_pp = pd.read_csv(args.promoter_peak)

    # subset shared observations (samples)
    common_obs = sorted(set(adata_rna.obs_names).intersection(set(adata_atac.obs_names)))
    if len(common_obs) == 0:
        raise SystemExit('No overlapping obs between RNA and ATAC')

    rna = adata_rna[common_obs, :].copy()
    atac = adata_atac[common_obs, :].copy()

    # Build candidate pairs from promoter-peak CSV (no DEG filtering)
    candidate_pairs = df_pp[['peak_id', 'gene_name']].drop_duplicates().copy()
    candidate_pairs = candidate_pairs[candidate_pairs['peak_id'].isin(atac.var_names) & candidate_pairs['gene_name'].isin(rna.var_names)].reset_index(drop=True)

    if candidate_pairs.empty:
        raise SystemExit('No candidate pairs after intersecting with var_names')

    print(f'Total candidate pairs: {len(candidate_pairs)}')

    # Prepare matrices
    peak_names = candidate_pairs['peak_id'].unique()
    gene_names = candidate_pairs['gene_name'].unique()

    print('Extracting matrices...')
    peak_mat = atac[:, peak_names].X
    gene_mat = rna[:, gene_names].X
    peak_mat = to_dense(peak_mat)
    gene_mat = to_dense(gene_mat)

    # map vectors
    peak_value_map = {peak: peak_mat[:, i] for i, peak in enumerate(peak_names)}
    gene_value_map = {gene: gene_mat[:, i] for i, gene in enumerate(gene_names)}

    # chunk and parallelize
    nchunks = max(1, args.chunks)
    chunks = [chunk.copy() for chunk in np.array_split(candidate_pairs, nchunks) if len(chunk) > 0]

    print(f'Processing in {len(chunks)} chunks using up to {min(args.chunks, len(chunks))} threads')
    all_results = []
    with ThreadPoolExecutor(max_workers=min(args.chunks, len(chunks))) as ex:
        for res in ex.map(lambda c: process_chunk(c, peak_value_map, gene_value_map, args.min_samples), chunks):
            all_results.extend(res)

    df_corr = pd.DataFrame(all_results)
    full_out = out_prefix.with_suffix('.peak_gene_corr.full.csv')
    sig_out = out_prefix.with_suffix('.peak_gene_corr.sig.csv')

    if df_corr.empty:
        print('No correlations computed. Writing empty files.')
        df_corr.to_csv(full_out, index=False)
        pd.DataFrame(columns=['peak_id', 'gene_name', 'correlation', 'pval', 'n_samples', 'FDR']).to_csv(sig_out, index=False)
        return

    print('Applying FDR correction...')
    df_corr['FDR'] = multipletests(df_corr['pval'].values, method='fdr_bh')[1]
    df_corr = df_corr.sort_values(['FDR', 'pval', 'correlation'], ascending=[True, True, False]).reset_index(drop=True)
    df_corr.to_csv(full_out, index=False)

    df_sig = df_corr[df_corr['FDR'] < 0.05].copy().reset_index(drop=True)
    df_sig.to_csv(sig_out, index=False)

    print('Wrote:', full_out)
    print('Wrote:', sig_out)


if __name__ == '__main__':
    main()
