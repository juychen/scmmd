import os
import sys
import subprocess
import multiprocessing
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed

# --------------------------------------------------------
# 日志配置
# --------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("/data1st1/junyi/output/sn0827/scenic_auc/aucell_processing.log")
    ]
)
logger = logging.getLogger(__name__)

# --------------------------------------------------------
# 处理单个 AUCell 任务
# --------------------------------------------------------
def run_aucell(ctx_file, expr_dir, ctx_dir, auc_dir, workers):
    """运行单个 AUCell 任务"""
    base_name = ctx_file.split('ctx_')[-1].split('.csv')[0]
    expr_file = os.path.join(expr_dir, f"{base_name}.loom")
    ctx_path = os.path.join(ctx_dir, ctx_file)
    auc_output = os.path.join(auc_dir, f"auc_{base_name}.tsv")
    flag_file = os.path.join(auc_dir, f"{base_name}.processing")

    # 跳过已处理文件
    if os.path.exists(auc_output):
        return f"Skip (exists): {base_name}"
    if os.path.exists(flag_file):
        return f"Skip (processing): {base_name}"

    # 检查表达矩阵
    if not os.path.exists(expr_file):
        return f"Missing expr: {expr_file}"

    # 写入 flag
    with open(flag_file, "w") as f:
        f.write("processing\n")

    cmd = [
        "pyscenic", "aucell",
        expr_file,
        ctx_path,
        "--output", auc_output,
        "--num_workers", str(workers)
    ]

    logger.info(f"▶ Running AUCell for: {base_name}")

    try:
        subprocess.run(cmd, check=True)
        logger.info(f"✅ Completed AUCell: {base_name}")
    except subprocess.CalledProcessError as e:
        logger.error(f"❌ Failed {base_name}: {e}")
    finally:
        if os.path.exists(flag_file):
            os.remove(flag_file)

    return base_name

# --------------------------------------------------------
# 主函数：并行执行
# --------------------------------------------------------
def parallel_aucell(expr_dir, ctx_dir, auc_dir, workers, parallel_jobs=4):

    os.makedirs(auc_dir, exist_ok=True)

    ctx_files = [f for f in os.listdir(ctx_dir) if f.startswith("ctx_") and f.endswith(".csv")]

    # 预筛选：只提交输出缺失且无 flag 的任务
    files_to_run = []
    for ctx_file in ctx_files:
        base_name = ctx_file.split('ctx_')[-1].split('.csv')[0]
        auc_output = os.path.join(auc_dir, f"auc_{base_name}.tsv")
        flag_file = os.path.join(auc_dir, f"{base_name}.processing")
        if os.path.exists(auc_output):
            print(f"Skip: {base_name} (output exists)")
        elif os.path.exists(flag_file):
            print(f"Skip: {base_name} (flag exists)")
        else:
            files_to_run.append(ctx_file)

    logger.info(f"🚀 Submitting {len(files_to_run)}/{len(ctx_files)} AUCell tasks...")
    logger.info(f"👉 Parallel jobs: {parallel_jobs}")
    logger.info(f"👉 Workers per job (--num_workers): {workers}")

    results = []

    with ThreadPoolExecutor(max_workers=parallel_jobs) as executor:
        future_map = {
            executor.submit(run_aucell, ctx_file, expr_dir, ctx_dir, auc_dir, workers): ctx_file
            for ctx_file in files_to_run
        }

        for future in as_completed(future_map):
            ctx_file = future_map[future]
            try:
                result = future.result()
                logger.info(f"Done: {result}")
                results.append(result)
            except Exception as e:
                logger.error(f"⚠️ Exception for {ctx_file}: {e}")

    logger.info("🎉 All AUCell tasks completed.")
    return results


# --------------------------------------------------------
# 运行入口
# --------------------------------------------------------
if __name__ == "__main__":
    EXPR_DIR = "/data1st1/junyi/output/sn0827/"
    CTX_DIR = "/data1st1/junyi/output/sn0827/scenic_ctx/"
    AUC_DIR = "/data1st1/junyi/output/sn0827/scenic_auc/"

    workers = min(32, max(1, multiprocessing.cpu_count() - 4))

    parallel_aucell(EXPR_DIR, CTX_DIR, AUC_DIR, workers, parallel_jobs=4)
