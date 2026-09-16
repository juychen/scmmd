#!/bin/bash
source /home/junyichen/anaconda3/etc/profile.d/conda.sh
conda activate allcools
export NUMEXPR_MAX_THREADS=40
set -euo pipefail
shopt -s nullglob

# 输入：L1subtypes/ 下每个 L1 细胞类型一个 .txt，每行一个 allc 路径
L1_DIR=/data2st1/junyi/methlyatlas/mCseq/CEMBA.mC.Metadata/L1subtypes
# 输出：<L1>.tsv.gz（每个细胞类型一个）、<Region>.tsv.gz（每个脑区一个）
OUT_DIR=/data2st1/junyi/merge_allcnew

# 每个 merge 任务内部用的线程数
CPU="${CPU:-32}"
# 同时并行运行的任务数。4 个任务 = 4x32 = 128 线程（本机 384 核）。
# 想更快可以调大，例如 PARALLEL=8 CPU=32（256 核）。
# 注意负载很不均：HPF_Glut 一个就占了 51367 个输入里的 25749 个，
# 所以 4 路并行下它单独仍是关键路径（约 30-40 小时），调大 PARALLEL 帮不到它。
PARALLEL="${PARALLEL:-4}"
# 每个任务一个日志文件，免得 4 路并行的输出互相穿插
LOG_DIR="$OUT_DIR/logs"
# DRY_RUN=1：只打印将要执行的命令，不写任何文件
DRY_RUN="${DRY_RUN:-0}"

mkdir -p "$OUT_DIR"

# chrom size 文件由你自己提供（两列：chr 前缀的染色体名 <TAB> 长度）。
# 脚本不生成它，只检查是否存在；想换个位置就改这一行。
CHROM_SIZE=/data2st1/junyi/single_allcs/mm10.main.chrom.sizes
if [[ ! -f "$CHROM_SIZE" ]]; then
    echo "ERROR: chrom size file not found: $CHROM_SIZE" >&2
    echo "       put your own two-column chrom size file there, or edit CHROM_SIZE in this script." >&2
    exit 1
fi

if [[ "$DRY_RUN" == 1 ]]; then
    FAILED=/dev/null
else
    mkdir -p "$LOG_DIR"
    FAILED="$LOG_DIR/failed.txt"
    : > "$FAILED"

    # 跳过判断只看输出文件是否存在，所以两个实例会把彼此还没建立占位文件的任务
    # 各跑一遍、互相覆盖。发现已有 allcools merge-allc 在跑就直接退出。
    if pgrep -f 'merge-allc' >/dev/null 2>&1; then
        echo "ERROR: 已有 allcools merge-allc 进程在运行，拒绝启动第二个实例：" >&2
        pgrep -af 'merge-allc' >&2
        echo "       请等它结束（或先 kill）再重跑本脚本。" >&2
        exit 1
    fi
fi

ts() { date '+%F %T'; }

# 输出文件存在就跳过 —— 包括 0 字节的占位文件。allcools 一开始就建占位文件、
# 最后才写 .tbi，所以占位文件说明该任务正在跑或上次被中断；两种情况脚本都不去
# 覆盖它，要重跑请手动删掉那个文件。返回 0 表示应当跳过。
should_skip() {   # $1=输出路径  $2=标签
    local out="$1" tag="$2"
    [[ -e "$out" ]] || return 1
    if [[ -f "$out.tbi" ]]; then
        echo "[skip] $tag (already exists)"
    else
        echo "[skip] $tag (占位文件：可能正在跑或上次中断，不覆盖；要重跑请手动删除 $out)" >&2
    fi
    return 0
}

# 任务池：最多同时跑 PARALLEL 个任务
RUNNING=0
wait_for_slot() {
    if (( RUNNING >= PARALLEL )); then
        wait -n 2>/dev/null || true
        RUNNING=$(( RUNNING - 1 ))
    fi
}

# run_merge <tag> <allc_paths_list> <output>
# 失败只记进 $FAILED，不让 set -e 把整轮跑杀掉：一个坏细胞类型不该终止其余 21 个任务。
run_merge() {
    local tag="$1" paths="$2" out="$3"
    if [[ "$DRY_RUN" == 1 ]]; then
        echo "$(ts) [dry-run] allcools merge-allc --cpu $CPU --chrom_size_path $CHROM_SIZE --allc_paths $paths --output_path $out"
        return 0
    fi
    if allcools merge-allc --cpu "$CPU" \
            --chrom_size_path "$CHROM_SIZE" \
            --allc_paths "$paths" \
            --output_path "$out" > "$LOG_DIR/$tag.log" 2>&1; then
        echo "$(ts) [done] $tag"
    else
        echo "$(ts) [fail] $tag (see $LOG_DIR/$tag.log)" >&2
        echo "$tag" >> "$FAILED"
    fi
}

# ---- Step 1: 每个 L1 细胞类型 merge 成一个文件（最多 PARALLEL 个并行）----
for list in "$L1_DIR"/*.txt; do
    l1=$(basename "$list" .txt)
    out="$OUT_DIR/$l1.tsv.gz"
    if should_skip "$out" "$l1"; then
        continue
    fi
    wait_for_slot
    echo "$(ts) [merge] $l1 (cpu=$CPU)"
    run_merge "$l1" "$list" "$out" &
    RUNNING=$(( RUNNING + 1 ))
done
wait
RUNNING=0

# ---- Step 2: 等 Step 1 完成后，每个脑区把所有细胞类型 merge 成一个文件 ----
for region in $(for f in "$L1_DIR"/*.txt; do basename "$f" .txt; done | cut -d_ -f1 | sort -u); do
    l1_lists=("$L1_DIR/${region}"_*.txt)
    expected=${#l1_lists[@]}

    if [[ "$DRY_RUN" == 1 ]]; then
        merge_list=$(mktemp -t merge_allc_dryrun.XXXXXX)
    else
        merge_list="$OUT_DIR/${region}_L1list.txt"
    fi
    : > "$merge_list"
    got=0
    for f in "$OUT_DIR/${region}"_*.tsv.gz; do
        # allcools 的批次中间文件（<out>.gzbatch_N.tmp.tsv.gz）不是细胞类型，别混进来
        [[ "$f" == *.gzbatch_* ]] && continue
        # 只收完整输出：没有 .tbi 说明还没写完，收进去会污染 region 结果
        [[ -f "$f.tbi" ]] || continue
        echo "$f" >> "$merge_list"
        got=$(( got + 1 ))
    done

    if (( got < expected )); then
        echo "[skip] region $region (可用 $got / 应有 $expected 个细胞类型，先补齐 Step 1)" >&2
        continue
    fi
    out="$OUT_DIR/$region.tsv.gz"
    if should_skip "$out" "$region"; then
        continue
    fi
    wait_for_slot
    echo "$(ts) [merge] region $region ($got cell types)"
    run_merge "$region" "$merge_list" "$out" &
    RUNNING=$(( RUNNING + 1 ))
done
wait

if [[ -s "$FAILED" ]]; then
    echo "=== $(wc -l < "$FAILED") 个任务失败 ===" >&2
    while IFS= read -r t; do
        echo "  - $t" >&2
    done < "$FAILED"
    echo "重跑本脚本即可：已完成的任务会被跳过，失败的任务会重试。" >&2
    exit 1
fi

echo "done. outputs in $OUT_DIR"
