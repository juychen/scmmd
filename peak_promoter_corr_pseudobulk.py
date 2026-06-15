#!/usr/bin/env python3
"""
Compute peak-promoter ATAC Pearson correlations across pseudobulk samples.

For each peak–promoter pair (from a 500 kb window), the script correlates the
ATAC signal of the distal peak with the ATAC signal of the linked promoter
across shared pseudobulk samples.

Usage:
  python peak_promoter_corr_pseudobulk.py \
    --adata-atac /path/to/adata_atac.h5ad \
    --promoter-peak /path/to/promoter_peak_500kb.csv \
    --out-prefix /path/to/output_prefix \
    --chunks 20

The promoter-peak CSV should contain at least:
    peak_id    — distal peak identifier (must be in adata_atac.var_names)
    promoter   — promoter peak identifier (must be in adata_atac.var_names)
Optionally columns gene_name and distance will be propagated to output.
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
    p = argparse.ArgumentParser(
        description="Peak–Promoter ATAC correlation on pseudobulk samples"
    )
    p.add_argument(
        '--adata-atac', required=True,
        help='h5ad file for pseudobulk ATAC (used for both peak and promoter values)'
    )
    p.add_argument(
        '--promoter-peak', required=True,
        help='CSV with columns [peak_id, promoter, …] from bedtools closest / window'
    )
    p.add_argument(
        '--promoter-col', default='promoter',
        help='Column name in --promoter-peak for the promoter peak ID (default: promoter)'
    )
    p.add_argument(
        '--out-prefix', required=True,
        help='Output prefix for result files'
    )
    p.add_argument(
        '--chunks', type=int, default=20,
        help='Number of chunks/threads (default 20)'
    )
    p.add_argument(
        '--min-samples', type=int, default=3,
        help='Minimum samples required to compute a correlation (default 3)'
    )
    return p.parse_args()


def to_dense(x):
    if sparse.issparse(x):
        return np.asarray(x.toarray())
    return np.asarray(x)


def process_chunk(chunk_df, peak_value_map, prom_value_map, min_samples=3):
    """Correlate each (peak, promoter) pair across samples."""
    results = []
    for row in chunk_df.itertuples(index=False):
        peak = row.peak_id
        prom = row.promoter
        pv = peak_value_map.get(peak)
        pr = prom_value_map.get(prom)
        if pv is None or pr is None:
            continue
        pv = np.asarray(pv).ravel()
        pr = np.asarray(pr).ravel()
        mask = np.isfinite(pv) & np.isfinite(pr)
        pv2 = pv[mask]
        pr2 = pr[mask]
        if pv2.size < min_samples:
            continue
        if np.all(pv2 == pv2[0]) or np.all(pr2 == pr2[0]):
            continue
        try:
            r, pval = pearsonr(pv2, pr2)
        except Exception:
            continue
        # Propagate gene_name and distance if available
        out = {
            'peak_id': peak,
            'promoter': prom,
            'correlation': float(r),
            'pval': float(pval),
            'n_samples': int(pv2.size),
        }
        if hasattr(row, 'gene_name') and not pd.isna(row.gene_name):
            out['gene_name'] = row.gene_name
        if hasattr(row, 'distance') and not pd.isna(row.distance):
            out['distance'] = row.distance
        results.append(out)
    return results


def main():
    args = parse_args()
    out_prefix = Path(args.out_prefix)
    prom_col = args.promoter_col

    print('Loading ATAC data...')
    adata_atac = sc.read_h5ad(args.adata_atac)

    print('Loading promoter–peak table...')
    df_pp = pd.read_csv(args.promoter_peak)

    # Validate required columns
    for col in ['peak_id', prom_col]:
        if col not in df_pp.columns:
            raise KeyError(f"Column '{col}' not found in {args.promoter_peak}")

    # Build candidate pairs
    keep_cols = ['peak_id', prom_col]
    extra_cols = [c for c in ['gene_name', 'distance'] if c in df_pp.columns]
    keep_cols += extra_cols

    candidate_pairs = df_pp[keep_cols].drop_duplicates().copy()
    candidate_pairs = candidate_pairs.rename(columns={prom_col: 'promoter'})

    # ---- Convert underscore IDs (chr1_1_100) to colon-dash (chr1:1-100) ----
    import bisect
    import re

    _region_re = re.compile(r'^(chr\w+):(\d+)-(\d+)$')

    def parse_region(s):
        """'chr1:1-100' -> ('chr1', 1, 100); invalid -> None"""
        m = _region_re.match(str(s).strip())
        if m:
            return (m.group(1), int(m.group(2)), int(m.group(3)))
        return None

    def underscore_to_region(s):
        """chr1_1_100 -> chr1:1-100"""
        parts = str(s).split('_')
        if len(parts) >= 3 and parts[0].startswith('chr'):
            return f"{parts[0]}:{parts[1]}-{parts[2]}"
        return str(s)

    # Build interval index from adata_atac.var_names: chrom -> [(start, end, var_name), ...] sorted by start
    chrom_intervals = {}
    for vn in adata_atac.var_names:
        parsed = parse_region(vn)
        if parsed is None:
            continue
        chrom, start, end = parsed
        chrom_intervals.setdefault(chrom, []).append((start, end, vn))

    for chrom in chrom_intervals:
        chrom_intervals[chrom].sort(key=lambda x: x[0])

    def find_overlap_var(chrom, start, end):
        """Return the var_name of first interval that overlaps (start, end), or None."""
        ivs = chrom_intervals.get(chrom)
        if not ivs:
            return None
        idx = bisect.bisect_left(ivs, (start, None, None), key=lambda x: x[0])
        if idx > 0 and ivs[idx - 1][1] > start:
            return ivs[idx - 1][2]
        while idx < len(ivs) and ivs[idx][0] <= end:
            if ivs[idx][1] > start:
                return ivs[idx][2]
            idx += 1
        return None

    # Convert and resolve: map CSV IDs → actual var_names via interval overlap
    candidate_pairs['promoter_region'] = candidate_pairs['promoter'].map(underscore_to_region)
    candidate_pairs['peak_region'] = candidate_pairs['peak_id'].map(underscore_to_region)

    resolved_peaks = []
    resolved_proms = []
    for _, row in candidate_pairs.iterrows():
        pp = parse_region(row['promoter_region'])
        pk = parse_region(row['peak_region'])
        prom_var = find_overlap_var(*pp) if pp else None
        peak_var = find_overlap_var(*pk) if pk else None
        resolved_peaks.append(peak_var)
        resolved_proms.append(prom_var)

    candidate_pairs['peak_id'] = resolved_peaks
    candidate_pairs['promoter'] = resolved_proms
    candidate_pairs = candidate_pairs.drop(columns=['promoter_region', 'peak_region'])
    candidate_pairs = candidate_pairs.dropna(subset=['peak_id', 'promoter']).reset_index(drop=True)

    if candidate_pairs.empty:
        raise SystemExit('No candidate pairs after intersecting with ATAC var_names')

    print(f'Total candidate pairs: {len(candidate_pairs)}')

    # Extract ATAC values for all unique peaks and promoters
    peak_ids = candidate_pairs['peak_id'].unique()
    prom_ids = candidate_pairs['promoter'].unique()
    all_ids = list(set(peak_ids).union(set(prom_ids)))

    print(f'Unique peak IDs: {len(peak_ids)}, unique promoter IDs: {len(prom_ids)}')

    print('Extracting ATAC matrix...')
    atac_mat = adata_atac[:, all_ids].X
    atac_mat = to_dense(atac_mat)

    # Build lookup maps
    peak_value_map = {}
    prom_value_map = {}
    for i, pid in enumerate(all_ids):
        vec = atac_mat[:, i]
        if pid in peak_ids:
            peak_value_map[pid] = vec
        if pid in prom_ids:
            prom_value_map[pid] = vec

    # Chunk and parallelise
    nchunks = max(1, args.chunks)
    chunks = [
        chunk.copy()
        for chunk in np.array_split(candidate_pairs, nchunks)
        if len(chunk) > 0
    ]

    print(f'Processing in {len(chunks)} chunks using up to {min(args.chunks, len(chunks))} threads')
    all_results = []
    with ThreadPoolExecutor(max_workers=min(args.chunks, len(chunks))) as ex:
        for res in ex.map(
            lambda c: process_chunk(c, peak_value_map, prom_value_map, args.min_samples),
            chunks,
        ):
            all_results.extend(res)

    df_corr = pd.DataFrame(all_results)

    # Output column order
    out_columns = ['peak_id', 'promoter']
    for c in ['gene_name', 'distance']:
        if c in df_corr.columns:
            out_columns.append(c)
    out_columns += ['correlation', 'pval', 'n_samples']

    full_out = out_prefix.with_suffix('.peak_promoter_corr.full.csv')
    sig_out = out_prefix.with_suffix('.peak_promoter_corr.sig.csv')

    if df_corr.empty:
        print('No correlations computed. Writing empty files.')
        pd.DataFrame(columns=out_columns + ['FDR']).to_csv(full_out, index=False)
        pd.DataFrame(columns=out_columns + ['FDR']).to_csv(sig_out, index=False)
        return

    print('Applying FDR correction...')
    df_corr['FDR'] = multipletests(df_corr['pval'].values, method='fdr_bh')[1]
    df_corr = df_corr.sort_values(
        ['FDR', 'pval', 'correlation'],
        ascending=[True, True, False],
    ).reset_index(drop=True)

    df_corr[out_columns + ['FDR']].to_csv(full_out, index=False)

    df_sig = df_corr[df_corr['FDR'] < 0.05].copy().reset_index(drop=True)
    df_sig[out_columns + ['FDR']].to_csv(sig_out, index=False)

    print('Wrote:', full_out)
    print('Wrote:', sig_out)


if __name__ == '__main__':
    main()
