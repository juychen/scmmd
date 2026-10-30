#!/usr/bin/env python3
"""PySCENIC CTX (cisTarget regulon pruning) for the region-level looms.

Same scheme as `ctx_sn0827.py`:
    pyscenic ctx adj_<region>.tsv <rankings.feather>
        --annotations_fname <motifs.tbl>
        --expression_mtx_fname <region>.loom
        --output ctx_<region>.csv --num_workers N

Output : <INPUT_DIR>/scenic_ctx/ctx_<region>.csv
Logs   : <INPUT_DIR>/scenic_ctx/ctx_run.log and logs/<region>.ctx.log
"""

import logging
import multiprocessing
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

# === PATHS ===
INPUT_DIR = "/data1st1/junyi/output/sn1029Region"
GRN_DIR = os.path.join(INPUT_DIR, "scenic_grn")
CTX_DIR = os.path.join(INPUT_DIR, "scenic_ctx")
LOG_DIR = os.path.join(CTX_DIR, "logs")
ANNOTATIONS = "/data1st2/yejun/scenic/mouse/motifs-v9-nr.mgi-m0.001-o0.0.tbl"
DATABASE = ("/data1st2/yejun/scenic/mouse/"
            "mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather")

# === RESOURCES ===
PARALLEL_JOBS = 4
NUM_WORKERS = min(32, max(1, multiprocessing.cpu_count() - 4))

os.makedirs(CTX_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(os.path.join(CTX_DIR, "ctx_run.log")),
    ],
)
logger = logging.getLogger(__name__)


def fix_grn_columns(grn_path):
    """Ensure the GRN file has the columns PySCENIC expects (TF, target, importance)."""
    try:
        with open(grn_path, "r") as fh:
            columns = fh.readline().strip().split("\t")
        target_columns = ["TF", "target", "importance"]
        if columns != target_columns:
            df = pd.read_csv(grn_path, sep="\t")
            df.columns = [c.lower() for c in df.columns]
            df = df.rename(columns={"tf": "TF"})
            df = df[["TF", "target", "importance"]]
            df.to_csv(grn_path, sep="\t", index=False)
            logger.info("✅ fixed GRN columns: %s", os.path.basename(grn_path))
        return True
    except Exception as exc:
        logger.error("⚠️ column fix failed for %s: %s", grn_path, exc)
        return False


def run_ctx(loom_file):
    base = os.path.splitext(loom_file)[0]
    expr_path = os.path.join(INPUT_DIR, loom_file)
    grn_path = os.path.join(GRN_DIR, f"adj_{base}.tsv")
    output_path = os.path.join(CTX_DIR, f"ctx_{base}.csv")
    flag_file = os.path.join(CTX_DIR, f"{base}.processing")
    job_log = os.path.join(LOG_DIR, f"{base}.ctx.log")

    if os.path.exists(flag_file):
        logger.info("⏩ skip (processing flag): %s", base)
        return f"skip-flag: {base}"
    if os.path.exists(output_path) and os.path.getsize(output_path) > 1024:
        logger.info("✅ skip (output exists): %s", base)
        return f"skip-exists: {base}"
    if not os.path.exists(grn_path):
        logger.warning("🚫 skip (GRN missing): %s", base)
        return f"missing-grn: {base}"

    with open(flag_file, "w") as fh:
        fh.write(f"Processing {loom_file}\n")

    if not fix_grn_columns(grn_path):
        os.remove(flag_file)
        return f"bad-grn: {base}"

    cmd = [
        "pyscenic", "ctx",
        grn_path,
        DATABASE,
        "--annotations_fname", ANNOTATIONS,
        "--expression_mtx_fname", expr_path,
        "--output", output_path,
        "--num_workers", str(NUM_WORKERS),
    ]

    t0 = time.time()
    logger.info("▶ CTX start: %s (workers=%d)", base, NUM_WORKERS)
    try:
        with open(job_log, "w") as fh:
            subprocess.run(cmd, check=True, stdout=fh, stderr=subprocess.STDOUT)
        size_mb = os.path.getsize(output_path) / 1e6
        logger.info("✅ CTX done : %s (%.1f min, %.1f MB)", base, (time.time() - t0) / 60, size_mb)
        return f"done: {base}"
    except subprocess.CalledProcessError as exc:
        logger.error("❌ CTX FAILED: %s (%s) -- see %s", base, exc, job_log)
        return f"failed: {base}"
    finally:
        if os.path.exists(flag_file):
            os.remove(flag_file)


def main():
    loom_files = sorted(f for f in os.listdir(INPUT_DIR) if f.endswith(".loom"))
    if not loom_files:
        logger.error("no .loom in %s", INPUT_DIR)
        return 1

    logger.info("=" * 70)
    logger.info("PySCENIC CTX -- %d looms", len(loom_files))
    logger.info("parallel jobs=%d, --num_workers=%d", PARALLEL_JOBS, NUM_WORKERS)
    logger.info("database: %s", DATABASE)
    logger.info("motifs  : %s", ANNOTATIONS)
    logger.info("=" * 70)

    with ThreadPoolExecutor(max_workers=PARALLEL_JOBS) as executor:
        futures = {executor.submit(run_ctx, f): f for f in loom_files}
        for fut in as_completed(futures):
            loom = futures[fut]
            try:
                logger.info("result: %s", fut.result())
            except Exception as exc:
                logger.error("⚠️ %s raised: %s", loom, exc)

    n_done = len([f for f in os.listdir(CTX_DIR) if f.startswith("ctx_") and f.endswith(".csv")])
    logger.info("🎉 CTX stage finished -- %d/%d ctx_*.csv present", n_done, len(loom_files))
    return 0


if __name__ == "__main__":
    sys.exit(main())
