#!/usr/bin/env python3
"""PySCENIC AUCell for the region-level looms.

Same scheme as `aucell_sn0827.py`:
    pyscenic aucell <region>.loom ctx_<region>.csv --output auc_<region>.tsv --num_workers N

Output : <INPUT_DIR>/scenic_auc/auc_<region>.tsv
Logs   : <INPUT_DIR>/scenic_auc/auc_run.log and logs/<region>.auc.log
"""

import logging
import multiprocessing
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

# === PATHS ===
INPUT_DIR = "/data1st1/junyi/output/sn1029Region"
CTX_DIR = os.path.join(INPUT_DIR, "scenic_ctx")
AUC_DIR = os.path.join(INPUT_DIR, "scenic_auc")
LOG_DIR = os.path.join(AUC_DIR, "logs")

# === RESOURCES ===
PARALLEL_JOBS = 4
NUM_WORKERS = min(32, max(1, multiprocessing.cpu_count() - 4))

os.makedirs(AUC_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(os.path.join(AUC_DIR, "auc_run.log")),
    ],
)
logger = logging.getLogger(__name__)


def run_aucell(loom_file):
    base = os.path.splitext(loom_file)[0]
    expr_path = os.path.join(INPUT_DIR, loom_file)
    ctx_path = os.path.join(CTX_DIR, f"ctx_{base}.csv")
    output_path = os.path.join(AUC_DIR, f"auc_{base}.tsv")
    flag_file = os.path.join(AUC_DIR, f"{base}.processing")
    job_log = os.path.join(LOG_DIR, f"{base}.auc.log")

    if os.path.exists(flag_file):
        logger.info("⏩ skip (processing flag): %s", base)
        return f"skip-flag: {base}"
    if os.path.exists(output_path) and os.path.getsize(output_path) > 1024:
        logger.info("✅ skip (output exists): %s", base)
        return f"skip-exists: {base}"
    if not os.path.exists(ctx_path):
        logger.warning("🚫 skip (CTX missing): %s", base)
        return f"missing-ctx: {base}"

    with open(flag_file, "w") as fh:
        fh.write(f"Processing {loom_file}\n")

    cmd = [
        "pyscenic", "aucell",
        expr_path,
        ctx_path,
        "--output", output_path,
        "--num_workers", str(NUM_WORKERS),
    ]

    t0 = time.time()
    logger.info("▶ AUCell start: %s (workers=%d)", base, NUM_WORKERS)
    try:
        with open(job_log, "w") as fh:
            subprocess.run(cmd, check=True, stdout=fh, stderr=subprocess.STDOUT)
        size_mb = os.path.getsize(output_path) / 1e6
        logger.info("✅ AUCell done : %s (%.1f min, %.1f MB)", base, (time.time() - t0) / 60, size_mb)
        return f"done: {base}"
    except subprocess.CalledProcessError as exc:
        logger.error("❌ AUCell FAILED: %s (%s) -- see %s", base, exc, job_log)
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
    logger.info("PySCENIC AUCell -- %d looms", len(loom_files))
    logger.info("parallel jobs=%d, --num_workers=%d", PARALLEL_JOBS, NUM_WORKERS)
    logger.info("=" * 70)

    with ThreadPoolExecutor(max_workers=PARALLEL_JOBS) as executor:
        futures = {executor.submit(run_aucell, f): f for f in loom_files}
        for fut in as_completed(futures):
            loom = futures[fut]
            try:
                logger.info("result: %s", fut.result())
            except Exception as exc:
                logger.error("⚠️ %s raised: %s", loom, exc)

    n_done = len([f for f in os.listdir(AUC_DIR) if f.startswith("auc_") and f.endswith(".tsv")])
    logger.info("🎉 AUCell stage finished -- %d/%d auc_*.tsv present", n_done, len(loom_files))
    return 0


if __name__ == "__main__":
    sys.exit(main())
