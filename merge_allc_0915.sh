#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate allcools
export NUMEXPR_MAX_THREADS=40
set -euo pipefail

# 输入：L1subtypes/ 下每个 L1 细胞类型一个 .txt，每行一个 allc 路径
L1_DIR=/data2st1/junyi/methlyatlas/mCseq/CEMBA.mC.Metadata/L1subtypes
# 输出：<L1>.tsv.gz（每个细胞类型一个）、<Region>.tsv.gz（每个脑区一个）
OUT_DIR=/data2st1/junyi/merge_allcnew
CPU=32

mkdir -p "$OUT_DIR"

# chrom size 文件由你自己提供（两列：chr 前缀的染色体名 <TAB> 长度）。
# 脚本不生成它，只检查是否存在；想换个位置就改这一行。
CHROM_SIZE=/data2st1/junyi/single_allcs/mm10.main.chrom.sizes
if [[ ! -f "$CHROM_SIZE" ]]; then
    echo "ERROR: chrom size file not found: $CHROM_SIZE" >&2
    echo "       put your own two-column chrom size file there, or edit CHROM_SIZE in this script." >&2
    exit 1
fi

# ---- Step 1: 每个 L1 细胞类型 merge 成一个文件 ----
for list in "$L1_DIR"/*.txt; do
    [[ -e "$list" ]] || continue
    l1=$(basename "$list" .txt)
    out="$OUT_DIR/$l1.tsv.gz"
    if [[ -f "$out" ]]; then
        echo "[skip] $l1 (already exists)"
        continue
    fi
    echo "[merge] $l1"
    allcools merge-allc --cpu "$CPU" \
        --chrom_size_path "$CHROM_SIZE" \
        --allc_paths "$list" \
        --output_path "$out"
done

# ---- Step 2: 等 Step 1 完成后，每个脑区把所有细胞类型 merge 成一个文件 ----
for region in $(for f in "$L1_DIR"/*.txt; do basename "$f" .txt; done | cut -d_ -f1 | sort -u); do
    merge_list="$OUT_DIR/${region}_L1list.txt"
    : > "$merge_list"
    for f in "$OUT_DIR/${region}"_*.tsv.gz; do
        if [[ -f "$f" ]]; then
            echo "$f" >> "$merge_list"
        fi
    done
    if [[ ! -s "$merge_list" ]]; then
        echo "[skip] region $region (no per-celltype output)"
        continue
    fi
    out="$OUT_DIR/$region.tsv.gz"
    if [[ -f "$out" ]]; then
        echo "[skip] region $region (already exists)"
        continue
    fi
    echo "[merge] region $region ($(wc -l < "$merge_list") cell types)"
    allcools merge-allc --cpu "$CPU" \
        --chrom_size_path "$CHROM_SIZE" \
        --allc_paths "$merge_list" \
        --output_path "$out"
done

echo "done. outputs in $OUT_DIR"
