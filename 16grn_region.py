#!/usr/bin/env python3
"""PySCENIC GRN (adjacency) for the region-level looms in sn1029Region.

One `pyscenic grn` job per loom (AMY / PFC / TH / HIP / STR / HY / MB), N jobs in
parallel, each with --num_workers workers -- the same scheme as
`grn_sn0827.py` (which produced sn0827L1/scenic_grn/adj_*.tsv).

Output : <INPUT_DIR>/scenic_grn/adj_<region>.tsv
Logs   : <INPUT_DIR>/scenic_grn/grn_run.log and logs/<region>.grn.log
Skip   : output exists (>1 KB) or a .processing flag is present
"""

import logging
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

# === PATHS ===
INPUT_DIR = "/data1st1/junyi/output/sn1029Region"
GRN_DIR = os.path.join(INPUT_DIR, "scenic_grn")
LOG_DIR = os.path.join(GRN_DIR, "logs")
TFS_FILE = "/data1st2/yejun/scenic/mouse/mm_mgi_tfs.txt"

# === RESOURCES ===
NUM_WORKERS = 32        # --num_workers per pyscenic grn job
PARALLEL_JOBS = 7       # loom files processed at the same time

os.makedirs(GRN_DIR, exist_ok=True)
os.makedirs(LOG_DIR, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(os.path.join(GRN_DIR, "grn_run.log")),
    ],
)
logger = logging.getLogger(__name__)


def run_grn(loom_file):
    """Run `pyscenic grn` on one loom, with a .processing flag file."""
    base = os.path.splitext(loom_file)[0]
    expr_path = os.path.join(INPUT_DIR, loom_file)
    output_path = os.path.join(GRN_DIR, f"adj_{base}.tsv")
    flag_file = os.path.join(GRN_DIR, f"{base}.processing")
    job_log = os.path.join(LOG_DIR, f"{base}.grn.log")

    if os.path.exists(flag_file):
        logger.info("⏩ skip (processing flag): %s", base)
        return f"skip-flag: {base}"
    if os.path.exists(output_path) and os.path.getsize(output_path) > 1024:
        logger.info("✅ skip (output exists): %s", base)
        return f"skip-exists: {base}"

    with open(flag_file, "w") as fh:
        fh.write(f"Processing {loom_file}\n")

    cmd = [
        "pyscenic", "grn",
        "--num_workers", str(NUM_WORKERS),
        "-o", output_path,
        expr_path,
        TFS_FILE,
    ]

    t0 = time.time()
    logger.info("▶ GRN start: %s (workers=%d)", base, NUM_WORKERS)
    try:
        with open(job_log, "w") as fh:
            subprocess.run(cmd, check=True, stdout=fh, stderr=subprocess.STDOUT)
        size_mb = os.path.getsize(output_path) / 1e6
        logger.info("✅ GRN done : %s (%.1f min, %.1f MB)", base, (time.time() - t0) / 60, size_mb)
        return f"done: {base}"
    except subprocess.CalledProcessError as exc:
        logger.error("❌ GRN FAILED: %s (%s) -- see %s", base, exc, job_log)
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
    logger.info("PySCENIC GRN -- %d looms in %s", len(loom_files), INPUT_DIR)
    logger.info("parallel jobs=%d, --num_workers=%d, TFs=%s",
                PARALLEL_JOBS, NUM_WORKERS, TFS_FILE)
    logger.info("looms: %s", ", ".join(os.path.splitext(f)[0] for f in loom_files))
    logger.info("=" * 70)

    with ThreadPoolExecutor(max_workers=PARALLEL_JOBS) as executor:
        futures = {executor.submit(run_grn, f): f for f in loom_files}
        for fut in as_completed(futures):
            loom = futures[fut]
            try:
                logger.info("result: %s", fut.result())
            except Exception as exc:
                logger.error("⚠️ %s raised: %s", loom, exc)

    n_done = len([f for f in os.listdir(GRN_DIR) if f.startswith("adj_") and f.endswith(".tsv")])
    logger.info("🎉 GRN stage finished -- %d/%d adj_*.tsv present", n_done, len(loom_files))
    return 0


if __name__ == "__main__":
    sys.exit(main())
