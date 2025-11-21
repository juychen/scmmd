import os
import sys
from operator import index
import pandas as pd
import numpy as np
import scanpy as sc
from typing import Optional, List, Dict, Literal
import glob

import subprocess
def run_r_script_deg_stacked_bars(file_path, sexes,  mode="manual",output_dir=None):
    if mode == "manual":
        script_path = "/data2st2/junyi/code/sn/example_deg_sex_direction_bar_function.R"
    elif mode == "auto":
        script_path = "/data2st2/junyi/code/sn/example_deg_sex_direction_bar_function_auto.R"
    else:
        raise ValueError("mode must be 'manual' or 'auto'")
    cmd = [
        "/home/junyichen/anaconda3/envs/r_env/bin/Rscript",
        script_path,
        file_path,
        ",".join(sexes)
    ]
    if output_dir:
        cmd.append(output_dir)
    try:
        subprocess.run(cmd, check=True)
        print("R script finished successfully.")
    except subprocess.CalledProcessError as e:
        print("Error running R script:", e)
        sys.exit(1)

def process_cell_counts(
    adata,
    region_order,
    region_col='region',
    celltype_col='celltype',
    gender_col='gender',
    status_col='status',
    neurotransmitter_col='Neurotransmitter',  # None
):
    """
    Output = [Summary rows] + [Per-region detail rows]
      • Summary: if a celltype spans multiple regions -> region='all' and counts summed;
                 if only one region -> keep that region.
      • Details: per-region rows for multi-region celltypes (only if neurotransmitter_col is not None).
    Row sort: by [region, celltype] only.
    If neurotransmitter_col is None, skip neurotransmitter logic and details.
    """
    # build gender_status
    adata.obs['gender_status'] = (
        adata.obs[gender_col].astype(str) + '_' + adata.obs[status_col].astype(str)
    )

    # crosstab: (region, celltype) x gender_status
    cross_tab = pd.crosstab(
        index=[adata.obs[region_col].astype('string'),
               adata.obs[celltype_col].astype('string')],
        columns=adata.obs['gender_status'].astype('string'),
        margins=False
    ).rename_axis(index=[region_col, celltype_col], columns=None).reset_index()

    # numeric columns
    value_cols = [c for c in cross_tab.columns if c not in [region_col, celltype_col]]
    cross_tab[value_cols] = cross_tab[value_cols].fillna(0).astype('int64')

    # prepare neurotransmitter mapping
    ct2nt, ct2nt4 = {}, {}
    if neurotransmitter_col is not None:
        nt_df = (adata.obs[[celltype_col, neurotransmitter_col]]
                 .astype({celltype_col: 'string', neurotransmitter_col: 'string'})
                 .dropna(subset=[celltype_col, neurotransmitter_col])
                 .drop_duplicates())
        dup_counts = nt_df.groupby(celltype_col)[neurotransmitter_col].nunique()
        if (dup_counts > 1).any():
            raise ValueError(
                "Each celltype must map to exactly one neurotransmitter; "
                f"violations for: {', '.join(map(str, dup_counts[dup_counts>1].index.tolist()))}"
            )
        ct2nt = dict(zip(nt_df[celltype_col], nt_df[neurotransmitter_col]))
        ct2nt4 = {k: (v[:4] if isinstance(v, str) else v) for k, v in ct2nt.items()}

    # regions covered by each celltype
    region_rank = {r: i for i, r in enumerate(region_order)}

    def _sorted_unique_regions(s):
        regs = list(dict.fromkeys(s))
        regs = [r for r in regs if pd.notna(r)]
        return sorted(regs, key=lambda x: region_rank.get(x, 10**9))

    ct2regions = cross_tab.groupby(celltype_col)[region_col].apply(_sorted_unique_regions)
    multi_cts = ct2regions.index[ct2regions.map(len) > 1]

    # ----- Part 1: Summary -----
    agg = cross_tab.groupby(celltype_col, as_index=False)[value_cols].sum()

    def _agg_region_for_ct(ct):
        regs = ct2regions.loc[ct]
        if neurotransmitter_col is None:
            return 'all' if len(regs) > 1 else (regs[0] if regs else 'all')
        nt = ct2nt.get(ct, None)
        if len(regs) > 1:
            return 'all'
        if len(regs) == 1 and nt == "NN":
            return 'all'
        return regs[0] if regs else 'all'

    agg[region_col] = agg[celltype_col].map(_agg_region_for_ct)

    # category order for region
    present_regions_agg = [r for r in region_order if r in set(agg[region_col])]
    if 'all' in set(agg[region_col]) and 'all' not in present_regions_agg:
        present_regions_agg = present_regions_agg + ['all']
    agg[region_col] = pd.Categorical(agg[region_col], categories=present_regions_agg, ordered=True)
    agg_sorted = agg.sort_values([region_col, celltype_col], kind='stable').reset_index(drop=True)

    # ----- Part 2: Details -----
    if neurotransmitter_col is None:
        detail_sorted = pd.DataFrame(columns=agg_sorted.columns)
    else:
        detail_cts = cross_tab[celltype_col].unique()
        nn_cts = [ct for ct, nt in ct2nt.items() if nt == "NN"]
        detail_cts = set(multi_cts).union(nn_cts)
        detail = cross_tab[cross_tab[celltype_col].isin(detail_cts)].copy()
        present_regions_detail = [r for r in region_order if r in set(detail[region_col])]
        detail[region_col] = pd.Categorical(detail[region_col], categories=present_regions_detail, ordered=True)
        detail_sorted = detail.sort_values([region_col, celltype_col], kind='stable').reset_index(drop=True)

    # ----- Concatenate (Summary first, then Details) -----
    out = pd.concat([agg_sorted, detail_sorted], ignore_index=True)

    # Reorder columns
    front = [region_col, celltype_col]
    rest = [c for c in out.columns if c not in front]
    out = out[front + rest]

    # Add neurotransmitter col if needed
    if neurotransmitter_col is not None:
        out[neurotransmitter_col] = out[celltype_col].map(ct2nt4)
        front = [region_col, neurotransmitter_col, celltype_col]
        rest = [c for c in out.columns if c not in front]
        out = out[front + rest]

    return out

def process_deg_counts(
    deg_path,
    celltype_col='celltype',
    gender_col='gender',
    status_col='status',
    region_col=None,
):
    """
    Count DEGs per (region?, celltype, gender) with columns ordered as:
    F_Down, F_Up, F_total, M_Down, M_Up, M_total (after index columns).
    """
    # read
    if deg_path.endswith('.csv'):
        df = pd.read_csv(deg_path, index_col=0)
    elif deg_path.endswith('.tsv') or deg_path.endswith('.txt'):
        df = pd.read_csv(deg_path, sep='\t', index_col=0)
    elif deg_path.endswith('.xlsx'):
        df = pd.read_excel(deg_path, index_col=0)
    else:
        raise ValueError("Unsupported file format; must be .csv, .tsv, .txt, or .xlsx")

    # required columns
    needed = [celltype_col, gender_col, status_col]
    if region_col is not None:
        needed.append(region_col)
    missing = [c for c in needed if c not in df.columns]
    if missing:
        try:
            df['Sex'] = 'M'
            df['Direction'] = df['status']
            df['Subclass'] = df['celltype.L2']
            df['Region'] = df['region']
        except KeyError:
            raise ValueError(f"Missing required columns in input: {missing}")

    # normalize relevant columns
    for c in needed:
        df[c] = df[c].astype(str).str.strip()

    # gender -> M/F (case-insensitive); keep other values unchanged
    g_lower = df[gender_col].str.lower()
    df[gender_col] = np.where(
        g_lower.isin(['male', 'm']), 'M',
        np.where(g_lower.isin(['female', 'f']), 'F', df[gender_col])
    )

    # status helper (case-insensitive)
    df['_status_norm'] = df[status_col].str.strip().str.casefold()

    # grouping keys
    if region_col is None:
        group_keys = [celltype_col, gender_col]
        index_cols = [celltype_col]
    else:
        group_keys = [region_col, celltype_col, gender_col]
        index_cols = [region_col, celltype_col]

    # aggregate
    grouped = (
        df.groupby(group_keys, dropna=False)
          .agg(
              total=(status_col, 'size'),
              Up=('_status_norm', lambda s: (s == 'up').sum()),
              Down=('_status_norm', lambda s: (s == 'down').sum()),
          )
          .reset_index()
    )

    # long -> wide
    long = grouped.melt(
        id_vars=group_keys,
        var_name='metric',
        value_name='count'
    )
    # keep a stable metric order for the pivot; we will re-order explicitly later
    long['metric'] = pd.Categorical(long['metric'],
                                    categories=['total', 'Up', 'Down'],
                                    ordered=True)

    wide = (
        long.pivot_table(
            index=index_cols,
            columns=[gender_col, 'metric'],
            values='count',
            fill_value=0,
            aggfunc='sum'
        )
        .sort_index(axis=1)  # temporary; final order set below
    )

    # flatten ('F','Down') -> 'F_Down'
    def _fmt(x):
        return 'NA' if pd.isna(x) else str(x)
    wide.columns = [f"{_fmt(g)}_{m}" for g, m in wide.columns.to_flat_index()]
    wide = wide.reset_index()

    # enforce final column order
    preferred = ['F_Down', 'F_Up', 'F_total', 'M_Down', 'M_Up', 'M_total']
    present_pref = [c for c in preferred if c in wide.columns]

    # index columns first, then preferred, then any remaining metric columns
    front = index_cols
    remaining = [c for c in wide.columns if c not in front + present_pref]
    wide = wide[front + present_pref + remaining]

    # filter: keep only rows with at least 1 DEG
    deg_cols = [c for c in ['F_Up', 'F_Down', 'M_Up', 'M_Down'] if c in wide.columns]
    wide = wide[wide[deg_cols].sum(axis=1) > 0]

    return wide


def split_celltype_region(
        df: pd.DataFrame,
        celltype_col: str = "celltype",
        region_col: str = "region",
        new_celltype_col: Optional[str] = None,
        Neurotransmitter_col: Optional[str] = None,
        inplace: bool = False,
) -> pd.DataFrame:
    """
    Split the `celltype_col` on the **first** space into two columns:
      - `region_col`: text before the first space
      - `new_celltype_col` (or `celltype_col` if None): text after the first space

    If `Neurotransmitter_col` is provided:
      - For cells with Neurotransmitter_col == 'NN', keep original celltype
      - For others, extract region but keep original celltype

    If a value has no space, both `region_col` and the (new) celltype receive the original value.
    Preserves NaNs (does not coerce them to the string "nan").

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame.
    celltype_col : str
        Source column to split (e.g., "celltype").
    region_col : str
        Name of the region column to create/overwrite.
    new_celltype_col : str | None
        If None, overwrite `celltype_col`; otherwise write the remainder to this new column.
    Neurotransmitter_col : str | None
        If provided, cells with this column != 'NN' will keep original celltype.
    inplace : bool
        If True, modify `df` in place; otherwise return a copy.

    Returns
    -------
    pd.DataFrame
        DataFrame with `region_col` inserted as the first column and the split celltype column
        (overwritten or newly created) added.
    """
    if celltype_col not in df.columns:
        raise KeyError(f"Column `{celltype_col}` not found in DataFrame.")

    if Neurotransmitter_col is not None and Neurotransmitter_col not in df.columns:
        raise KeyError(f"Column `{Neurotransmitter_col}` not found in DataFrame.")

    out = df if inplace else df.copy()

    # Original values (keep NaNs as NA, not strings)
    s = out[celltype_col].astype("string").str.strip()

    # Split only at the first whitespace
    parts = s.str.extract(r"^(?P<_region>\S+)\s+(?P<_rest>.+)$")

    # Fallback for rows without a space: use original value for both
    region_series = parts["_region"].fillna(s)
    rest_series = parts["_rest"].fillna(s)

    # Write/overwrite the region column and move it to the first position
    out[region_col] = region_series
    cols = [region_col] + [c for c in out.columns if c != region_col]
    out = out.loc[:, cols]

    # Handle Neurotransmitter_col logic
    if Neurotransmitter_col is not None:
        # For cells with Neurotransmitter_col != 'NN', keep original celltype
        nn_mask = out[Neurotransmitter_col] == 'NN'

        # Initialize target column with original values
        target_celltype_col = new_celltype_col or celltype_col
        out[target_celltype_col] = out[celltype_col]

        # Only update non-NN cells with the split values
        out.loc[nn_mask, target_celltype_col] = rest_series[nn_mask]
    else:
        # Original behavior: update all cells
        target_celltype_col = new_celltype_col or celltype_col
        out[target_celltype_col] = rest_series

    return out


def process_go_terms(
    go_file_path,
    sheet_name, # sheet name in the Excel file
    ontology_type,
    source_key='source',       # column name in the CSV file (e.g., 'GO')
    celltype_col='celltype',
    gender_col='sex',
    direction_col='direction',
    term_col='term_name',
    region_col=None,               # set to a column name to split by region; None to ignore
):
    """
    Count unique GO terms per (region?, celltype, sex_dir), where
    sex_dir = f"{ontology_type}_{<sex>}_{<direction>}".
    Returns a wide table with one column per sex_dir.
    """
    # read the sheet
    if go_file_path.endswith('.csv'):
        df = pd.read_csv(go_file_path)
    elif go_file_path.endswith('.xlsx'):
        df = pd.read_excel(go_file_path, sheet_name=sheet_name)
    else:
        raise ValueError("Unsupported file format; must be .csv or .xlsx")

    df = df[df[source_key]==f'GO:{ontology_type}']

    # validate columns
    needed = [celltype_col, gender_col, direction_col, term_col]
    if region_col is not None:
        needed.append(region_col)
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    # clean relevant columns (strings)
    for c in needed:
        df[c] = df[c].astype(str).str.strip()

    # gender -> M/F (case-insensitive); keep other values unchanged
    g_lower = df[gender_col].str.lower()
    df[gender_col] = np.where(
        g_lower.isin(['male', 'm']), 'M',
        np.where(g_lower.isin(['female', 'f']), 'F', df[gender_col])
    )

    # build sex_dir: "<ontology_type>_<sex>_<direction>"
    df['sex_dir'] = (
        str(ontology_type) + '_' +
        df[gender_col].astype(str).str.strip() + '_' +
        df[direction_col].astype(str).str.strip()
    )

    # choose index (with or without region)
    index_cols = [celltype_col] if region_col is None else [region_col, celltype_col]

    # pivot: count unique term names per index × sex_dir
    out = (
        df.pivot_table(
            index=index_cols,
            columns='sex_dir',
            values=term_col,
            aggfunc=lambda s: s.nunique(),
            fill_value=0
        )
        .sort_index(axis=1)
        .reset_index()
    )

    return out


def calculate_spatial_distribution(
        adata: sc.AnnData,
        region_key: str,
        subregion_key: str,
        celltype_key: str,
        mode: Literal['celltype_distribution', 'subregion_composition'] = 'celltype_distribution',
        group_key: Optional[str] = None,
        region_order: Optional[List[str]] = None,
        top_n: int = 1
):
    """
    Calculate spatial distribution patterns with region-based sorting

    Parameters:
        adata: AnnData object containing spatial data
        region_key: Column name for macro-regions
        subregion_key: Column name for sub-regions
        celltype_key: Column name for cell types
        mode: Analysis mode ('celltype_distribution' or 'subregion_composition')
        group_key: Optional grouping column
        region_order: Optional region ordering
        top_n: Number of top results to return

    Returns:
        DataFrame with analysis results sorted by region then celltype/subregion
    """
    # Validate input columns
    for key in [region_key, subregion_key, celltype_key]:
        if key not in adata.obs.columns:
            raise ValueError(f"Column '{key}' not found in adata.obs")

    # Prepare data with observed=True to avoid FutureWarning
    obs_df = adata.obs.copy()
    for col in obs_df.select_dtypes(include=['category']).columns:
        obs_df[col] = obs_df[col].astype(str)

    # Core calculation function
    def calculate_top_distribution(df, by_col: str, target_col: str):
        """Calculate percentage distribution between categories"""
        dist_df = df.groupby([by_col, target_col], observed=True).size().unstack(fill_value=0)
        dist_pct = dist_df.div(dist_df.sum(axis=1), axis=0) * 100

        if mode == 'celltype_distribution':
            result = pd.DataFrame({
                'DominantSubregion': dist_pct.idxmax(axis=1),
                'Percentage': dist_pct.max(axis=1).round(2)
            })
        else:
            result = pd.DataFrame()
            for idx in dist_pct.index:
                top_items = dist_pct.loc[idx].nlargest(top_n)
                for rank, (val, pct) in enumerate(zip(top_items.index, top_items), 1):
                    result.loc[idx, f'Top{rank}_CellType'] = val
                    result.loc[idx, f'Top{rank}_Percentage'] = round(pct, 2)
        return result, dist_pct

    # Run analysis based on mode
    try:
        if mode == 'celltype_distribution':
            main_result, _ = calculate_top_distribution(obs_df, celltype_key, subregion_key)
            main_result.index.name = 'CellType'
            main_result = main_result.reset_index()

            # Add region information
            if region_key:
                subregion_to_region = obs_df[[subregion_key, region_key]].drop_duplicates() \
                    .set_index(subregion_key)[region_key]

                # Handle categorical columns
                if isinstance(subregion_to_region, pd.Categorical):
                    subregion_to_region = subregion_to_region.astype(str)

                main_result['Region'] = main_result['DominantSubregion'].map(subregion_to_region)

                # Handle multi-region cell types
                celltype_regions = obs_df.groupby(celltype_key, observed=True)[region_key].unique()
                multi_region_mask = main_result['CellType'].isin(
                    celltype_regions[celltype_regions.apply(len) > 1].index
                )

                # Ensure 'Region' column can accept 'all' value
                if pd.api.types.is_categorical_dtype(main_result['Region']):
                    main_result['Region'] = main_result['Region'].astype(str)

                main_result.loc[multi_region_mask, 'Region'] = 'all'

        else:  # subregion_composition mode
            main_result, _ = calculate_top_distribution(obs_df, subregion_key, celltype_key)
            main_result.index.name = 'Subregion'
            main_result = main_result.reset_index()

            if region_key:
                subregion_to_region = obs_df[[subregion_key, region_key]].drop_duplicates() \
                    .set_index(subregion_key)[region_key]

                # Handle categorical columns
                if isinstance(subregion_to_region, pd.Categorical):
                    subregion_to_region = subregion_to_region.astype(str)

                main_result['Region'] = main_result['Subregion'].map(subregion_to_region)

    except KeyError as e:
        raise ValueError(f"Missing required data in adata.obs: {str(e)}")

    # Handle grouping if specified
    if group_key:
        if group_key not in adata.obs.columns:
            raise ValueError(f"Group column '{group_key}' not found")

        group_results = {}
        for group_name, group_df in obs_df.groupby(group_key, observed=True):
            group_result, _ = calculate_top_distribution(group_df,
                                                         celltype_key if mode == 'celltype_distribution' else subregion_key,
                                                         subregion_key if mode == 'celltype_distribution' else celltype_key)
            group_results[group_name] = group_result

        # Merge group results
        for group_name in group_results:
            if mode == 'celltype_distribution':
                main_result[f'{group_key}_{group_name}_Subregion'] = main_result['CellType'].map(
                    group_results[group_name]['DominantSubregion'])
                main_result[f'{group_key}_{group_name}_Pct'] = main_result['CellType'].map(
                    group_results[group_name]['Percentage'])
            else:
                for rank in range(1, top_n + 1):
                    main_result[f'{group_key}_{group_name}_Top{rank}_CellType'] = main_result['Subregion'].map(
                        group_results[group_name].get(f'Top{rank}_CellType', None))
                    main_result[f'{group_key}_{group_name}_Top{rank}_Pct'] = main_result['Subregion'].map(
                        group_results[group_name].get(f'Top{rank}_Percentage', None))

    # Apply region ordering if specified
    if region_order and region_key and 'Region' in main_result.columns:
        def get_region_order(region):
            if region == 'all':
                return len(region_order)
            try:
                return region_order.index(region)
            except ValueError:
                return len(region_order) + 1

        main_result['sort_key'] = main_result['Region'].apply(get_region_order)

        # Sort by region first, then by celltype/subregion alphabetically
        if mode == 'celltype_distribution':
            main_result = main_result.sort_values(['sort_key', 'CellType'])
        else:
            main_result = main_result.sort_values(['sort_key', 'Subregion'])

        main_result = main_result.drop(columns=['sort_key'])
    else:
        # Default alphabetical sorting
        if mode == 'celltype_distribution':
            main_result = main_result.sort_values('CellType')
        else:
            main_result = main_result.sort_values('Subregion')

    # Clean up column order
    if mode == 'celltype_distribution':
        column_order = ['Region', 'CellType', 'DominantSubregion', 'Percentage']
        if group_key:
            for group_name in obs_df[group_key].unique():
                column_order.extend([
                    f'{group_key}_{group_name}_Subregion',
                    f'{group_key}_{group_name}_Pct'
                ])
    else:
        column_order = ['Region', 'Subregion']
        for rank in range(1, top_n + 1):
            column_order.extend([f'Top{rank}_CellType', f'Top{rank}_Percentage'])
        if group_key:
            for group_name in obs_df[group_key].unique():
                for rank in range(1, top_n + 1):
                    column_order.extend([
                        f'{group_key}_{group_name}_Top{rank}_CellType',
                        f'{group_key}_{group_name}_Top{rank}_Pct'
                    ])

    # Replace NaN values before returning
    result_df = main_result[column_order].reset_index(drop=True)
    result_df = result_df.fillna({
        **{col: 0 for col in result_df.select_dtypes(include=['number']).columns},
        **{col: 'NA' for col in result_df.select_dtypes(exclude=['number']).columns}
    })

    return result_df



if __name__ == "__main__":
    class DummyAdata:
        def __init__(self, obs: pd.DataFrame):
            self.obs = obs

    # Load the dataset snRNA
    # adata_sn_path_dict = {
    #     '/data7/mark/STG/dataset/snRNA/merge_SCH_new/SUS_4v4_500_1000gene/adata_SUS_4v4_downsampled_ratio.h5ad':'SUS',
    #     '/data7/mark/STG/dataset/snRNA/merge_SCH_new/SUS_3v3_500_1000gene/adata_SUS_3v3_downsampled_ratio.h5ad':'SUS',
    #     "/data7/mark/STG/dataset/snRNA/merge_SCH_new/RES_4v3_500_1000gene/adata_RES_4v3_downsampled_ratio.h5ad":'RES',
    #     '/data7/mark/STG/dataset/snRNA/merge_SCH_new/RES_3v3_500_1000gene/adata_RES_3v3_downsampled_ratio.h5ad':'RES',
    #     '/data7/mark/STG/dataset/snRNA/merge_SCH_new/RES_3v3_500_1000gene_beirui/adata_RES_3v3_beirui_downsampled_ratio.h5ad':'RES',
    #     '/data7/mark/STG/dataset/snRNA/merge_SCH_new/CSDS_4v3_500_1000gene/adata_CSDS_4v3_downsampled_ratio.h5ad':'CSDS',
    #     '/data7/mark/STG/dataset/snRNA/merge_SCH_new/CSDS_3v3_500_1000gene/adata_CSDS_3v3_downsampled_ratio.h5ad':'CSDS',
    #     "/data7/mark/STG/dataset/snRNA/merge_SCH_new/CSDS_3v3_500_1000gene_beirui/adata_CSDS_3v3_beirui_downsampled_ratio.h5ad":'CSDS',
    #     '/data7/mark/STG/dataset/snRNA/merge_SCH_new/CSRES_4v3_500_1000gene/adata_CSRES_4v3_downsampled_ratio.h5ad':'CSRES',
    #     '/data7/mark/STG/dataset/snRNA/merge_SCH_new/CSRES_3v3_500_1000gene/adata_CSRES_3v3_downsampled_ratio.h5ad':'CSRES',
    #     '/data7/mark/STG/dataset/snRNA/merge_SCH_new/CSRES_3v3_500_1000gene_beirui/adata_CSRES_3v3_beirui_downsampled_ratio.h5ad':'CSRES'
    # }
    # adata_sn_path_dict = {
    #     '/data7/mark/STG/dataset/snRNA/merge_SCH_new/SUS_4v4_500_1000gene/SUS_4V4.h5ad':'SUS',
    #     '/data7/mark/STG/dataset/snRNA/merge_SCH_new/SUS_3v3_500_1000gene/SUS_3V3.h5ad':'SUS',
    #     "/data7/mark/STG/dataset/snRNA/merge_SCH_new/RES_4v3_500_1000gene/RES_4V3.h5ad":'RES',
    #     '/data7/mark/STG/dataset/snRNA/merge_SCH_new/RES_3v3_500_1000gene/RES_3V3.h5ad':'RES',
    #     '/data7/mark/STG/dataset/snRNA/merge_SCH_new/RES_3v3_500_1000gene_beirui/RES_3V3_beirui.h5ad':'RES',
    #     '/data7/mark/STG/dataset/snRNA/merge_SCH_new/CSDS_4v3_500_1000gene/CSDS_4V3.h5ad':'CSDS',
    #     '/data7/mark/STG/dataset/snRNA/merge_SCH_new/CSDS_3v3_500_1000gene/CSDS_3V3.h5ad':'CSDS',
    #     "/data7/mark/STG/dataset/snRNA/merge_SCH_new/CSDS_3v3_500_1000gene_beirui/CSDS_3V3_beirui.h5ad":'CSDS',
    #     '/data7/mark/STG/dataset/snRNA/merge_SCH_new/CSRES_4v3_500_1000gene/CSRES_4V3.h5ad':'CSRES',
    #     '/data7/mark/STG/dataset/snRNA/merge_SCH_new/CSRES_3v3_500_1000gene/CSRES_3V3.h5ad':'CSRES',
    #     '/data7/mark/STG/dataset/snRNA/merge_SCH_new/CSRES_3v3_500_1000gene_beirui/CSRES_3V3_beirui.h5ad':'CSRES'
    # }
    adata_sn_path_dict = {
        # '/data2st2/junyi/output/stg1028/CURES_4VN/CURES_4VN.h5ad':'RES',
        # '/data2st2/junyi/output/stg1028/CURES_3VN/CURES_3VN.h5ad':'RES',
        # '/data2st2/junyi/output/stg1028/CURES_3VB/CURES_3VB.h5ad':'RES',
        # '/data2st2/junyi/output/stg1028/CSRES_4VN/CSRES_4VN.h5ad':'CSRES',
        # '/data2st2/junyi/output/stg1028/CSRES_3VN/CSRES_3VN.h5ad':'CSRES',
        # '/data2st2/junyi/output/stg1028/CSRES_3VB/CSRES_3VB.h5ad':'CSRES',
        '/data2st1/junyi/output/atac1112/dar/celltype.L2/MASTNG_dar_annotated.csv':'SUS',
        # '/data2st2/junyi/output/stg1028/CSDS_3VN/CSDS_3VN.h5ad':'CSDS',
        # '/data2st2/junyi/output/stg1028/CSDS_3VB/CSDS_3VB.h5ad':'CSDS',
        #'/data2st2/junyi/output/stg1028/CUMS_4VN/CUMS_4VN.h5ad':'SUS'
        # '/data2st2/junyi/output/stg1028/CUMS_3VN/CUMS_3VN.h5ad':'SUS'
    }

    region_order = ['AMY', 'HPF', 'PFC']
    NT = ['Glut', 'GABA', 'Hist', 'Chol', 'Dopa', 'Sero']
    adj_method_list = [''] # 'fdr','bonferroni'
    memento_region_name = 'company_ngeneson_sample'
    log2fc_cutoff_list = ['annoatated']
    csv_suffix = "" # "" for original csv; "filtered" for filtered csv
    deg_prefix_list = ['MASTNG_dar']
    glob_string = "MASTNG*annotated.csv"

    df_Class = pd.read_excel(
        '/data2st2/junyi/output/name_form_new.xlsx')
    subclass_to_class_map = df_Class.drop_duplicates(subset=['Subclass']).set_index('Subclass')[
        'Class_for_degs'].to_dict()

    suffix_list = [
        f'{prefix}_{method}_{cutoff}'
        for prefix in deg_prefix_list
        for method in adj_method_list
        for cutoff in log2fc_cutoff_list
    ]

    for adata_sn_path, disease in adata_sn_path_dict.items():

        # === process snRNA ===
        #adata_sn = sc.read_h5ad(adata_sn_path)
        meta_atac = pd.read_csv(f'/data2st1/junyi/output/atac1112/ATACSC_3REGION_ALL_L2annoated.csv', index_col=0)
        region_key = 'Region'
        cell_type_key = 'celltype'
        group_key = 'Condition'
        gender_key = 'sex'
        sample_key = 'sample'
        disease = disease
        neurotransmitter_key = 'Neurotransmitter'
        meta_atac[gender_key] = "M"
        # remove cell that celltype.L1 =='OPC' and Neurotransmitter!='NN'
        meta_atac = meta_atac[~((meta_atac['celltype.L1_ct'] == 'OPC') & (meta_atac['Neurotransmitter_celltype'] != 'NN'))]
        if cell_type_key not in meta_atac:
            meta_atac[cell_type_key] = meta_atac['celltype.L2']
        if neurotransmitter_key not in meta_atac:
            meta_atac[neurotransmitter_key] = meta_atac['Neurotransmitter_celltype'].str[:4]
        if gender_key not in meta_atac:
            meta_atac.rename(columns={'gender': gender_key}, inplace=True)
        if 'CON' not in meta_atac[group_key].unique():
            meta_atac[group_key] = meta_atac[group_key].replace({'MW': 'CON', 'MC': disease})

        meta_atac['celltype.L1'] = meta_atac['celltype.L1'].astype(str)
        meta_atac.loc[meta_atac['Neurotransmitter'] != 'NN', 'celltype.L1'] = \
            meta_atac.loc[meta_atac['Neurotransmitter'] != 'NN', 'Neurotransmitter']
        meta_atac.loc[
            meta_atac['celltype.L1'] == 'Non-neuron', 'celltype.L1'] = 'Vascular'
        meta_atac.loc[
            meta_atac['celltype.L2'].isin(['COP', 'NFOL', 'MFOL']), 'celltype.L1'] = 'Oligo'
        meta_atac.loc[meta_atac['sample'].str.startswith('WT'), 'status'] = 'CON'

        # if 'beirui' in adata_sn_path or 'CSDS' in adata_sn_path or 'CSRES' in adata_sn_path:
        #     adata_sn = adata_sn[meta_atac[gender_key] == 'M'].copy()
        #adata_sn.write_h5ad(adata_sn_path)
        meta_atac[region_key] = meta_atac[region_key].str.replace('HIP', 'HPF', regex=False)
        valid_regions = [r for r in region_order if r in meta_atac[region_key].unique()]
        print(f"Valid regions: {valid_regions}")

        output_dir = os.path.dirname(adata_sn_path)
        dir_list = glob.glob(f"{output_dir}/{glob_string}")
        for dir in dir_list:
            dir_path = f'{output_dir}/{dir}'
            deg_candidates = [dir]  # base CSV path
            # if not os.path.isdir(dir_path):
            #     continue
            # # If csv_suffix is provided, add another DEG path
            # if csv_suffix:
            #     deg_candidates.append(f'{dir_path}/{dir}_{csv_suffix}.csv')
            
            for deg_path in deg_candidates:
                # Define output subfolder depending on which file is used
                # if deg_path.endswith(f'_{csv_suffix}.csv'):
                #     current_output_dir = f'{dir_path}/{dir}_{csv_suffix}'
                # else:
                current_output_dir = f'{dir}dir'
                os.makedirs(current_output_dir, exist_ok=True)

        #         # Check matching suffix
        #         for suffix in suffix_list:
        #             if dir.endswith(suffix) and os.path.exists(deg_path):

                for dir in dir_list:
                        adata_sn = DummyAdata(meta_atac)
                        valid_regions = [r for r in region_order if r in meta_atac[region_key].unique()]
                        sn_counts = process_cell_counts(
                            adata_sn,
                            region_order=valid_regions,
                            region_col=region_key,
                            celltype_col=cell_type_key,
                            gender_col=gender_key,
                            status_col=group_key,
                        )
                        sn_counts = sn_counts[sn_counts[region_key] != 'all']
                        sn_counts.rename(columns={region_key: 'Region', cell_type_key: 'Subclass'}, inplace=True)

                        # === cell type DEG number for snRNA ===
                        sn_deg_counts = process_deg_counts(deg_path,celltype_col='Subclass',gender_col='Sex',status_col='Direction',region_col='Region')
                        deg_cols = ["F_Down", "F_Up", "F_total", "M_Down", "M_Up", "M_total"]
                        deg_cols = [col for col in deg_cols if col in sn_deg_counts.columns]
                        # merge
                        sn_counts['Subclass'] = sn_counts['Subclass'].str.replace(" ", "_")
                        sn_counts['Subclass'] = sn_counts['Subclass'].str.replace("/", "-")
                        sn_merged = sn_counts.drop(columns=deg_cols, errors="ignore").merge(
                            sn_deg_counts[["Region", "Subclass"] + deg_cols],
                            on=["Region", "Subclass"],
                            how="left"
                        )
                        if deg_cols:
                            sn_merged[deg_cols] = sn_merged[deg_cols].fillna(0).astype(int)
                        if 'Class' not in sn_merged.columns:
                            sn_merged['Class'] = sn_merged['Subclass'].map(subclass_to_class_map)
                            cols = ['Class'] + [col for col in sn_merged.columns if col != 'Class']
                            sn_merged = sn_merged[cols]
                        sn_merged.to_csv(f'{current_output_dir}/3_sn_cell_DEG_number.csv', index=False)

                        # === statistics for snRNA ===
                        # sexes = meta_atac[gender_key].unique()
                        if any(x in current_output_dir for x in ['beirui', 'CSDS', 'CSRES','3VB']):
                            sexes = ['M']
                        else:
                            sexes = ['M']

                        sn_merged = sn_merged[sn_merged.isna().sum(axis=1) <= 10].reset_index(drop=True)
                        id_vars = ['Class', 'Subclass', 'Region', 'Neurotransmitter']
                        records = []
                        for _, row in sn_merged.iterrows():
                            for sex in sexes:
                                n_cells = row[f'{sex}_CON'] + row[f'{sex}_{disease}']
                                for direction in ['Down', 'Up']:
                                    n_degs = row[f'{sex}_{direction}']
                                    records.append({
                                        'Class': row['Class'],
                                        'Subclass': row['Subclass'],
                                        'Region': row['Region'],
                                        'Neurotransmitter': row['Neurotransmitter'],
                                        'Sex': sex,
                                        'Direction': direction,
                                        'n_cells': n_cells,
                                        'DEG_total': n_degs
                                    })

                        long_df = pd.DataFrame(records)
                        long_df['n_celltypes'] = 1
                        long_df.to_excel(f'{current_output_dir}/3_sn_cell_DEG_number_long_format.xlsx', index=False)

                        # plot bar figures
                        try:
                            run_r_script_deg_stacked_bars(
                                f'{current_output_dir}/3_sn_cell_DEG_number_long_format.xlsx',
                                sexes,
                                'manual',
                                f'{current_output_dir}/Figure2_deg_number_barplot',
                            )
                        except Exception as e:
                            print(f"Error in manual barplot for {current_output_dir}: {e}")

                        try:
                            run_r_script_deg_stacked_bars(
                                f'{current_output_dir}/3_sn_cell_DEG_number_long_format.xlsx',
                                sexes,
                                'auto',
                                f'{current_output_dir}/Figure2_deg_number_barplot_auto',
                            )
                        except Exception as e:
                            print(f"Error in auto barplot for {current_output_dir}: {e}")

                        # === deg statistics summary ===
                        if deg_path.endswith('.csv'):
                            deg_list = pd.read_csv(deg_path)
                        else:
                            deg_list = pd.read_excel(deg_path)
                        try:
                            deg_list['Sex'] = 'M'
                            deg_list['Direction'] = deg_list['status']
                            deg_list['Subclass'] = deg_list['celltype.L2']
                            deg_list['Region'] = deg_list['region']
                        except Exception as e:
                            raise ValueError(f"Missing required columns in input: {e}")

                        if 'Region subclass' not in deg_list.columns:
                            deg_list['Region subclass'] = deg_list.apply(
                                lambda row: f"{row['Region']} {row['Subclass']}" if row['Region'] not in row['Subclass'] else
                                row['Subclass'],
                                axis=1
                            )
                        if 'Class' not in deg_list.columns:
                            deg_list['Class'] = deg_list['Subclass'].map(subclass_to_class_map)

                        deg = deg_list.copy()
                        for col in ['Gene_name', 'Ensemble', 'Subclass', 'Sex', 'Direction']:
                            deg[col] = deg[col].astype(str).str.strip()


                        # 1. Preserve original statistics (without N/NN separation)
                        def create_original_stats(df, celltype_col='Region subclass', gene_col='Gene', ensemble_col='Ensemble',
                                                  sex_col='Sex',
                                                  direction_col='Direction'):
                            # Total cell type frequency (without gender)
                            ct_no_gender = (
                                df.groupby([gene_col, ensemble_col])[celltype_col]
                                .nunique()
                                .rename('ct_freq_total_no_gender')
                            )

                            # Cell type frequency including gender
                            df['celltype_sex'] = df[celltype_col] + '_' + df[sex_col]
                            ct_with_gender = (
                                df.groupby([gene_col, ensemble_col])['celltype_sex']
                                .nunique()
                                .rename('ct_freq_total')
                            )

                            # Group by gender and status
                            combo_grp = df.groupby([gene_col, ensemble_col, sex_col, direction_col])[celltype_col]
                            combo_freq = combo_grp.nunique().rename('freq')
                            combo_cts = combo_grp.apply(lambda s: ','.join(sorted(set(s)))).rename('cts')

                            # Convert to wide format
                            combo = pd.concat([combo_freq, combo_cts], axis=1).reset_index()
                            combo['col_freq'] = combo[sex_col] + '_' + combo[direction_col] + '_freq'
                            combo['col_cts'] = combo[sex_col] + '_' + combo[direction_col] + '_cts'

                            wide_freq = combo.pivot_table(
                                index=[gene_col, ensemble_col],
                                columns='col_freq',
                                values='freq',
                                fill_value=0,
                                aggfunc='first'
                            )

                            wide_cts = combo.pivot_table(
                                index=[gene_col, ensemble_col],
                                columns='col_cts',
                                values='cts',
                                fill_value='',
                                aggfunc='first'
                            )

                            return ct_no_gender, ct_with_gender, wide_freq, wide_cts


                        # 2. Create grouped statistics (separated by N and NN)
                        def create_grouped_stats(df_subset, suffix, celltype_col='Region subclass', gene_col='Gene',
                                                 ensemble_col='Ensemble',
                                                 sex_col='Sex', direction_col='Direction'):
                            # Cell type frequency (including gender)
                            ct_freq = (
                                df_subset.groupby([gene_col, ensemble_col])['celltype_sex']
                                .nunique()
                                .rename(f'ct_freq_total_{suffix}')
                            )

                            # Group statistics (gender + status)
                            combo_grp = df_subset.groupby([gene_col, ensemble_col, sex_col, direction_col])[celltype_col]
                            combo_freq = combo_grp.nunique().rename(f'freq_{suffix}')
                            combo_cts = combo_grp.apply(lambda s: ','.join(sorted(set(s)))).rename(f'cts_{suffix}')

                            # Convert to wide format
                            combo = pd.concat([combo_freq, combo_cts], axis=1).reset_index()
                            combo[f'col_freq_{suffix}'] = combo[sex_col] + '_' + combo[direction_col] + f'_freq_{suffix}'
                            combo[f'col_cts_{suffix}'] = combo[sex_col] + '_' + combo[direction_col] + f'_cts_{suffix}'

                            wide_freq = combo.pivot_table(
                                index=[gene_col, ensemble_col],
                                columns=f'col_freq_{suffix}',
                                values=f'freq_{suffix}',
                                fill_value=0,
                                aggfunc='first'
                            )

                            wide_cts = combo.pivot_table(
                                index=[gene_col, ensemble_col],
                                columns=f'col_cts_{suffix}',
                                values=f'cts_{suffix}',
                                fill_value='',
                                aggfunc='first'
                            )

                            return ct_freq, wide_freq, wide_cts


                        # Generate original statistics
                        orig_stats = create_original_stats(deg, celltype_col='Region subclass')

                        # Generate N/NN grouped statistics
                        stats_N = create_grouped_stats(deg[deg['Neurotransmitter'] != 'NN'], 'N',
                                                       celltype_col='Region subclass')
                        stats_NN = create_grouped_stats(deg[deg['Neurotransmitter'] == 'NN'], 'NN',
                                                        celltype_col='Region subclass')

                        # Combine all results
                        summary = pd.concat([
                            *orig_stats,  # Original statistics
                            *stats_N,  # N group statistics
                            *stats_NN  # NN group statistics
                        ], axis=1).reset_index()

                        # Sort and save
                        summary = summary.sort_values('ct_freq_total', ascending=False)
                        summary["Gene_ranking"] = range(1, len(summary) + 1)
                        cols = ["Gene_ranking"] + [c for c in summary.columns if c != "Gene_ranking"]
                        summary = summary[cols]
                        summary.to_csv(f"{current_output_dir}/1-2 all sn_DEG_list_summary with NN regions.csv",index=False)


                        def summarize_deg_specific_stats_by_range(df, column_name, thresholds):
                            """
                            Summarizes row counts and proportions for specific numerical ranges.

                            Parameters:
                            - df (pd.DataFrame): The DataFrame containing the data.
                            - column_name (str): The name of the column to perform the analysis on.
                            - thresholds (list): A sorted list of numerical values defining the ranges.
                                                 E.g., [0, 1, 5, 20] defines ranges [0, 1), [1, 5), [5, 20), and [20, ∞).

                            Returns:
                            - pd.DataFrame: A DataFrame with the count and proportion for each range.
                            """
                            if column_name not in df.columns:
                                print(f"Error: Column '{column_name}' not found in the DataFrame.")
                                return pd.DataFrame()

                            total_rows = len(df)
                            results = []

                            # Ensure thresholds are sorted
                            thresholds.sort()

                            # Iterate through the thresholds to define ranges
                            for i in range(len(thresholds)):
                                start_val = thresholds[i]

                                # Handle all intervals except the last one
                                if i < len(thresholds) - 1:
                                    end_val = thresholds[i + 1]
                                    # Count includes the start point but not the end point (>= start and < end)
                                    count = len(df[(df[column_name] > start_val) & (df[column_name] <= end_val)])
                                    label = f"({start_val}, {end_val}]"
                                # Handle the last interval (>= last value)
                                else:
                                    count = len(df[df[column_name] >= start_val])
                                    label = f"[{start_val}, ∞)"

                                proportion = count / total_rows if total_rows > 0 else 0
                                results.append({
                                    "Metric": label,
                                    "Count": count,
                                    "Total_Count": total_rows,
                                    "Proportion": proportion
                                })

                            return pd.DataFrame(results)

                        excel_file_name = f"{current_output_dir}/1-3 DEG_frequency_with_region_without_sex.xlsx"
                        with pd.ExcelWriter(excel_file_name) as writer:
                            # Summarize DEG-specific statistics by defined ranges
                            results_df = summarize_deg_specific_stats_by_range(summary, 'ct_freq_total_no_gender',
                                                                               [0, 1, 5, 20])
                            # Save to a specific sheet named after the key
                            results_df.to_excel(writer, index=False)


                        # === add class frequency to DEG summary ===
                        deg_list['Gene'] = deg_list['Gene_name'].str.replace('"', '')

                        # Group by gene, Sex, Class, and Direction, then count occurrences
                        direction_freq = (deg_list.groupby(['Gene', 'Sex', 'Class', 'Direction'])
                                          .size()
                                          .reset_index(name='count'))
                        direction_pivot = (direction_freq.pivot_table(
                            index='Gene',
                            columns=['Sex', 'Class', 'Direction'],
                            values='count',
                            fill_value=0
                        ))

                        direction_pivot.columns = [f"{sex}_{cls}_{dir}" for sex, cls, dir in direction_pivot.columns]
                        # Calculate total counts per Class (without separating by sex or direction)
                        class_total = (deg_list.groupby(['Gene', 'Class']).size().reset_index(name='count')
                                       .pivot_table(index='Gene', columns='Class', values='count', fill_value=0)
                                       )
                        # Calculate total counts per Sex and Class combination (without direction separation)
                        sex_class_total = (deg_list.groupby(['Gene', 'Sex', 'Class']).size().reset_index(name='count')
                                           .pivot_table(index='Gene', columns=['Sex', 'Class'], values='count', fill_value=0)
                                           )
                        sex_class_total.columns = [f"{sex}_{cls}_Total" for sex, cls in sex_class_total.columns]
                        # Merge all statistical results
                        final_freq = (direction_pivot.merge(sex_class_total, on='Gene', how='outer')
                                      .merge(class_total, on='Gene', how='outer', suffixes=('', '_Total'))
                                      .reset_index()
                                      )

                        # Define the desired column order based on classes
                        classes = ["Glut", "GABA", "Hist", "Chol", "Dopa", "Sero",
                                   "Astrocyte", "Microglia", "Oligo", "OPC", "Vascular", "Epen", "Immune"
                                   ]
                        ordered_cols = ['Gene']
                        for cls in classes:
                            if cls in final_freq.columns:
                                ordered_cols.append(cls)
                            for sex in sexes:
                                sex_cls_total_col = f"{sex}_{cls}_Total"
                                if sex_cls_total_col in final_freq.columns:
                                    ordered_cols.append(sex_cls_total_col)
                                for direction in ['Up', 'Down']:
                                    direction_col = f"{sex}_{cls}_{direction}"
                                    if direction_col in final_freq.columns:
                                        ordered_cols.append(direction_col)
                        final_freq = final_freq[ordered_cols]

                        # Merge with deg_summary data
                        deg_summary_final = summary.merge(final_freq, on='Gene', how='left')
                        deg_summary_final.to_csv(f"{current_output_dir}/1-2 all sn_DEG_list_summary with NN regions.csv", index=False)


                        # # === add go term summary ===
                        # duplicate_handlings = ["filter","Down","Up","nlogp_keep_highest","nlogp_compare"]
                        # for dup in duplicate_handlings:
                        #     go_path = f'{dir_path}/{dup}/2_go_merged_with_predictions.csv'
                        #     go_path_unique = f'{dir_path}/{dup}/2_go_merged_with_predictions_unique.csv'
                        #     if os.path.exists(go_path_unique):
                        #         go_path = go_path_unique
                        #
                        #     bp_summary = process_go_terms(go_path, ontology_type='BP',sheet_name='ALL',celltype_col='Subclass',
                        #                                   gender_col='Sex',direction_col='Direction',term_col='term_name',
                        #                                   region_col='Region'
                        #     )
                        #     mf_summary = process_go_terms(go_path, ontology_type='MF',sheet_name='ALL',celltype_col='Subclass',
                        #                                   gender_col='Sex',direction_col='Direction',term_col='term_name',
                        #                                   region_col='Region'
                        #     )
                        #     cc_summary = process_go_terms(go_path, ontology_type='CC',sheet_name='ALL',celltype_col='Subclass',
                        #                                   gender_col='Sex',direction_col='Direction',term_col='term_name',
                        #                                   region_col='Region'
                        #     )
                        #
                        #     final_df = (
                        #         sn_merged.merge(bp_summary, on=['Region','Subclass'], how='left')
                        #         .merge(mf_summary, on=['Region','Subclass'], how='left')
                        #         .merge(cc_summary, on=['Region','Subclass'], how='left')
                        #     )
                        #
                        #     # Calculate totals
                        #     for ont in ['BP', 'MF', 'CC']:
                        #         for sex in sexes:
                        #             up_col = f'{ont}_{sex}_Up'
                        #             down_col = f'{ont}_{sex}_Down'
                        #             tot_col = f'{ont}_{sex}_Total'
                        #             if up_col not in final_df.columns:   final_df[up_col] = 0
                        #             if down_col not in final_df.columns:   final_df[down_col] = 0
                        #             final_df[tot_col] = final_df[up_col].fillna(0).astype(float) + \
                        #                                 final_df[down_col].fillna(0).astype(float)
                        #
                        #     ordered_go_cols = []
                        #     for ont in ['BP', 'MF', 'CC']:
                        #         for sex in sexes:
                        #             for met in ['Down', 'Up', 'Total']:
                        #                 col = f'{ont}_{sex}_{met}'
                        #                 if col in final_df.columns:
                        #                     ordered_go_cols.append(col)
                        #
                        #     if ordered_go_cols:
                        #         final_df[ordered_go_cols] = (
                        #             final_df[ordered_go_cols]
                        #             .apply(pd.to_numeric, errors='coerce')  # in case some are object
                        #             .fillna(0)
                        #             .astype('int64')
                        #         )
                        #
                        #     other_cols = [c for c in final_df.columns if c not in set(ordered_go_cols)]
                        #     final_df = final_df[other_cols + ordered_go_cols]
                        #     cols = ['Class'] + [col for col in final_df.columns if col != 'Class']
                        #     final_df = final_df[cols]
                        #     # Save final DataFrame
                        #     final_df.to_csv(f"{dir_path}/3_sn_cell_DEG_Go_{dup}_number.csv", index=False)



    #
    # # ====== NN GO enrichment analysis for modules and DEG list ======
    # module_go_path = '/data7/mark/STG/dataset/ST_HD_Ex_out/function_analysis/sn_ST_cell_DEG_number/NN_new/2-2 NN_BP_MF_GO enrichment analysis.xlsx'
    # go_deg_path = '/data7/mark/STG/dataset/ST_HD_Ex_out/function_analysis/sn_ST_cell_DEG_number/NN_new/2-1 NN GO_orderT_p0.01_nlogP_noreplication_filtered with regions.csv'
    # deg_list_path = '/data7/mark/STG/dataset/ST_HD_Ex_out/function_analysis/sn_ST_cell_DEG_number/NN_new/1-1 NN DEG list MAST_Wilcoxon_fdr_logfc01_filtered with regions_new.csv'
    # # Process GO terms for modules
    # module_df = (pd.read_excel(module_go_path, sheet_name='NN_BP_MF_sorted',usecols=['term_name', 'module_manual']))
    #
    # # Process GO terms for DEG list
    # go_df = pd.read_csv(go_deg_path, index_col=0)
    # go_df = go_df[go_df['Neurotransmitter']=='NN'] # filter for non-neurotransmitter GO terms
    # go_df['term_name'] = go_df['source'] + '_' + go_df['term_name']
    # go_df_merged = go_df.merge(module_df, on='term_name', how='inner') # common go term
    #
    # go_expanded = (
    #     go_df_merged
    #     .assign(gene=lambda d: d['intersection'].str.split(','))
    #     .explode('gene')
    #     .assign(gene=lambda d: d['gene'].str.strip())  # clean whitespace
    # )
    # gene_sample_summary = (
    #     go_expanded
    #     .groupby(['gene', 'sample'], as_index=False)
    #     .agg(
    #         modules=('module_manual', lambda x: ';'.join(sorted(set(x)))),
    #         go_terms=('term_name', lambda x: ';'.join(sorted(set(x))))
    #     )
    # )
    #
    # deg_list = pd.read_csv(deg_list_path, index_col=0)
    # deg_list = deg_list[deg_list['Neurotransmitter'] == 'NN']  # filter for non-neurotransmitter DEG lists
    # deg_list = (
    #     deg_list
    #     .assign(
    #         gene=lambda d: d['gene_name'].str.replace(r'^"|"$', '', regex=True)
    #     )
    # )
    # deg_list['sample'] = deg_list['celltype'] + '_' + deg_list['gender']
    # meta_cols = ['gene', 'sample','gene_name', 'ensemble', 'avg_log2FC', 'p_val_adj', 'Neurotransmitter', 'region', 'gender', 'status',
    #    'Neurotransmitter_celltype', 'region_old']
    # deg_meta = deg_list[meta_cols].drop_duplicates(subset=['gene', 'sample'])
    #
    # gene_sample_summary = (
    #     gene_sample_summary
    #     .merge(deg_meta, on=['gene', 'sample'], how='left')
    # )
    # gene_sample_summary.rename(columns={'modules': 'Module', 'go_terms': 'GO'}, inplace=True)
    # gene_sample_summary.to_csv(os.path.join(output_dir_path,"NN_new", "GO_terms_connected_to_DEG_list.csv"), index=False)
    #
    #
    #
    # # ===== calculate GO term for FU Qinghui =====
    # module_go_path = '/data7/mark/STG/dataset/ST_HD_Ex_out/function_analysis/sn_ST_cell_DEG_number/N/2-2 GO enrichment analysis 20250804.xlsx'
    # go_deg_path = '/data7/mark/STG/dataset/ST_HD_Ex_out/function_analysis/sn_ST_cell_DEG_number/N/2-1 neuron GO list 20250729.xlsx'
    #
    # # Process GO terms for modules
    # module_df_highly = (pd.read_excel(module_go_path, sheet_name='N BP_MF_highly_ enriched V2',usecols=['Go', 'Module']))
    # module_df_low = (pd.read_excel(module_go_path, sheet_name='N BP_MF_low_enriched_V2',usecols=['Go', 'Module']))
    # module_df_CC = (pd.read_excel(module_go_path, sheet_name='N CC', usecols=['Go', 'module']))
    #
    #
    # # Process GO terms for DEG list
    # go_df = pd.read_excel(go_deg_path, sheet_name='ALL')
    # # go_df = go_df[go_df['Neurotransmitter']=='NN'] # filter for non-neurotransmitter GO terms
    # go_df['Go'] = go_df['source'] + '_' + go_df['term_name']
    #
    #
    # def pick_module(df, out_col):
    #     for c in ('module_manual', 'Module', 'module'):
    #         if c in df.columns:
    #             return df[['Go', c]].rename(columns={c: out_col})
    #     raise KeyError(f"No module column found in {out_col}")
    # mod_high = pick_module(module_df_highly, 'Module_high')
    # mod_low = pick_module(module_df_low, 'Module_low')
    # mod_cc = pick_module(module_df_CC, 'Module_CC')
    #
    # go_df_merged = (
    #     go_df.merge(mod_high, on='Go', how='left', validate='m:1')
    #     .merge(mod_low, on='Go', how='left', validate='m:1')
    #     .merge(mod_cc, on='Go', how='left', validate='m:1')
    # )
    # go_df_merged.to_excel(go_deg_path.replace('.xlsx', '_merged.xlsx'), index=False)
    #
    #

