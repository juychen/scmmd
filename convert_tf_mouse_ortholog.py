"""
convert_tf_mouse_ortholog.py

将 TF 名称转换为正式的老鼠基因名。
1. 先在老鼠基因中查
2. 如果不在老鼠中，到人基因中查
3. 如果在人基因中找到，通过 ortholog 转换为老鼠基因

Usage:
    python convert_tf_mouse_ortholog.py --tf-list-file <path> --output <path>
    python convert_tf_mouse_ortholog.py --tf-list <TF1,TF2,...> --output <path>
"""

import argparse
import pandas as pd
from pybiomart import Dataset


def convert_tf_to_mouse_ortholog(tf_list, dataset_mouse, dataset_human):
    """
    将 TF 名称转换为正式的老鼠基因名。

    Parameters:
        tf_list: list of TF names (uppercase)
        dataset_mouse: pybiomart Dataset for mouse (mmusculus_gene_ensembl)
        dataset_human: pybiomart Dataset for human (hsapiens_gene_ensembl)

    Returns:
        mapping_dict: {original_TF: converted_mouse_gene_name}
        summary: dict with counts
    """
    tf_set = set([tf.upper().strip() for tf in tf_list])

    # === Step 1: 下载所有老鼠基因名 ===
    print("Step 1: 下载老鼠基因名...")
    df_mouse = dataset_mouse.query(attributes=['external_gene_name', 'ensembl_gene_id'])
    mouse_gene_col = [c for c in df_mouse.columns if 'name' in c.lower()][0]
    mouse_genes = set(df_mouse[mouse_gene_col].dropna().str.upper().str.strip().unique())
    print(f"  老鼠基因总数: {len(mouse_genes)}")

    mouse_found = tf_set & mouse_genes
    mapping = {tf: tf for tf in mouse_found}
    unmatched = tf_set - mouse_found
    print(f"  老鼠中找到: {len(mouse_found)} 个")
    print(f"  未匹配: {len(unmatched)} 个")

    if not unmatched:
        return mapping, {'mouse_found': len(mouse_found), 'human_found': 0,
                         'human_to_mouse': 0, 'unmatched': 0}

    # === Step 2: 下载所有人类基因名 ===
    print("\nStep 2: 下载人类基因名...")
    df_human = dataset_human.query(attributes=['external_gene_name', 'ensembl_gene_id'])
    human_gene_col = [c for c in df_human.columns if 'name' in c.lower()][0]
    human_genes = set(df_human[human_gene_col].dropna().str.upper().str.strip().unique())
    print(f"  人类基因总数: {len(human_genes)}")

    human_matched = unmatched & human_genes
    still_unmatched = unmatched - human_genes
    print(f"  人类中找到: {len(human_matched)} 个")
    print(f"  仍未知: {len(still_unmatched)} 个")

    if not human_matched:
        for tf in still_unmatched:
            mapping[tf] = tf
        return mapping, {'mouse_found': len(mouse_found), 'human_found': len(human_matched),
                         'human_to_mouse': 0, 'unmatched': len(still_unmatched)}

    # === Step 3: 下载老鼠中所有人类同源基因映射 ===
    print(f"\nStep 3: 下载老鼠-人类同源基因映射表...")
    df_orth = dataset_mouse.query(
        attributes=['external_gene_name', 'hsapiens_homolog_associated_gene_name']
    )
    mouse_name_col = [c for c in df_orth.columns if 'gene' in c.lower() and 'name' in c.lower()][0]
    human_col = [c for c in df_orth.columns if 'human' in c.lower()][0]

    df_orth = df_orth.dropna(subset=[mouse_name_col, human_col])
    df_orth[mouse_name_col] = df_orth[mouse_name_col].astype(str).str.upper().str.strip()
    df_orth[human_col] = df_orth[human_col].astype(str).str.upper().str.strip()
    df_orth = df_orth[(df_orth[mouse_name_col] != 'NAN') & (df_orth[human_col] != 'NAN')]

    human2mouse = {}
    for _, row in df_orth.iterrows():
        hg = row[human_col]
        mg = row[mouse_name_col]
        if hg not in human2mouse:
            human2mouse[hg] = mg

    print(f"  人类→老鼠 同源映射数: {len(human2mouse)}")

    human_to_mouse = 0
    for tf in human_matched:
        if tf in human2mouse:
            mapping[tf] = human2mouse[tf]
            human_to_mouse += 1
        else:
            mapping[tf] = tf

    for tf in still_unmatched:
        mapping[tf] = tf

    print(f"\n=== 转换汇总 ===")
    print(f"老鼠基因名可直接使用: {len(mouse_found)}")
    print(f"人类基因转老鼠同源: {human_to_mouse}")
    print(f"人类基因但无同源: {len(human_matched) - human_to_mouse}")
    print(f"未找到(保留原名): {len(still_unmatched)}")

    return mapping, {
        'mouse_found': len(mouse_found),
        'human_found': len(human_matched),
        'human_to_mouse': human_to_mouse,
        'unmatched': len(still_unmatched)
    }


def main():
    parser = argparse.ArgumentParser(
        description="将 TF 名称转换为正式的老鼠基因名 (含人→鼠同源映射)"
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--tf-list-file', type=str,
                       help='包含 TF 列表的文件 (一列, 无表头)')
    group.add_argument('--tf-list', type=str,
                       help='逗号分隔的 TF 列表, 如 "AR,ATF4,BACH1"')
    parser.add_argument('--output', type=str, required=True,
                        help='输出 CSV 路径 (TF_original, TF_official 两列)')

    args = parser.parse_args()

    if args.tf_list_file:
        tf_list = pd.read_csv(args.tf_list_file, header=None).iloc[:, 0].dropna().tolist()
    else:
        tf_list = [t.strip() for t in args.tf_list.split(',') if t.strip()]

    print(f"加载了 {len(tf_list)} 个 TF")

    dataset_mouse = Dataset(name='mmusculus_gene_ensembl',
                            host='http://www.ensembl.org')
    dataset_human = Dataset(name='hsapiens_gene_ensembl',
                            host='http://www.ensembl.org')

    mapping, summary = convert_tf_to_mouse_ortholog(tf_list, dataset_mouse, dataset_human)

    df_out = pd.DataFrame({
        'TF_original': sorted(set([t.upper().strip() for t in tf_list]))
    })
    df_out['TF_official'] = df_out['TF_original'].map(mapping)
    df_out.to_csv(args.output, index=False)
    print(f"\n结果已保存至: {args.output}")


if __name__ == '__main__':
    main()
