#!/usr/bin/env python3
"""PySCENIC AUCell for the region-level looms.

Same scheme as `aucell_sn0827.py`:
    pyscenic aucell <region>.loom ctx_<region>.csv --output auc_<region>.tsv --num_workers N

Pass --ctx-dir/--output-dir to score a different regulon set (e.g. the v10nr
ctx run) without touching the default v9 results.

Output : <output-dir>/auc_<region>.tsv   (default <input-dir>/scenic_auc)
Logs   : <output-dir>/auc_run.log and logs/<region>.auc.log
"""

import argparse
import logging
import multiprocessing
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

# === DEFAULT PATHS ===
INPUT_DIR = "/data1st1/junyi/output/sn1029Region"
DEFAULT_PARALLEL = 4


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="PySCENIC AUCell on region-level looms")
    p.add_argument("--input-dir", default=INPUT_DIR, help="directory with the .loom files")
    p.add_argument("--ctx-dir", default=None, help="directory with ctx_*.csv (default <input-dir>/scenic_ctx)")
    p.add_argument("--output-dir", default=None, help="auc output directory (default <input-dir>/scenic_auc)")
    p.add_argument("--parallel", type=int, default=DEFAULT_PARALLEL, help="looms processed at once")
    p.add_argument("--workers", type=int, default=None, help="--num_workers (default: ncpu-4, max 32)")
    p.add_argument("--regions", default=None, help="comma separated subset, e.g. AMY,STR")
    return p.parse_args(argv)


def run_aucell(loom_file, cfg, logger):
    base = os.path.splitext(loom_file)[0]
    expr_path = os.path.join(cfg["input_dir"], loom_file)
    ctx_path = os.path.join(cfg["ctx_dir"], f"ctx_{base}.csv")
    output_path = os.path.join(cfg["output_dir"], f"auc_{base}.tsv")
    flag_file = os.path.join(cfg["output_dir"], f"{base}.processing")
    job_log = os.path.join(cfg["log_dir"], f"{base}.auc.log")

    if os.path.exists(flag_file):
        logger.info("⏩ skip (processing flag): %s", base)
        return f"skip-flag: {base}"
    if os.path.exists(output_path) and os.path.getsize(output_path) > 1024:
        logger.info("✅ skip (output exists): %s", base)
        return f"skip-exists: {base}"
    if not os.path.exists(ctx_path) or os.path.getsize(ctx_path) < 1024:
        logger.warning("🚫 skip (CTX missing/empty): %s", base)
        return f"missing-ctx: {base}"

    with open(flag_file, "w") as fh:
        fh.write(f"Processing {loom_file}\n")

    cmd = ["pyscenic", "aucell", expr_path, ctx_path,
           "--output", output_path,
           "--num_workers", str(cfg["workers"])]

    t0 = time.time()
    logger.info("▶ AUCell start: %s (workers=%d)", base, cfg["workers"])
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


def main(argv=None):
    args = parse_args(argv)

    input_dir = args.input_dir
    ctx_dir = args.ctx_dir or os.path.join(input_dir, "scenic_ctx")
    output_dir = args.output_dir or os.path.join(input_dir, "scenic_auc")
    log_dir = os.path.join(output_dir, "logs")
    workers = args.workers or min(32, max(1, multiprocessing.cpu_count() - 4))

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stdout),
                  logging.FileHandler(os.path.join(output_dir, "auc_run.log"))],
    )
    logger = logging.getLogger(__name__)

    loom_files = sorted(f for f in os.listdir(input_dir) if f.endswith(".loom"))
    if args.regions:
        wanted = {r.strip() for r in args.regions.split(",") if r.strip()}
        loom_files = [f for f in loom_files if os.path.splitext(f)[0] in wanted]
    if not loom_files:
        logger.error("no .loom in %s", input_dir)
        return 1

    cfg = dict(input_dir=input_dir, ctx_dir=ctx_dir, output_dir=output_dir,
               log_dir=log_dir, workers=workers)

    logger.info("=" * 70)
    logger.info("PySCENIC AUCell -- %d looms", len(loom_files))
    logger.info("parallel jobs=%d, --num_workers=%d", args.parallel, workers)
    logger.info("ctx dir   : %s", ctx_dir)
    logger.info("output    : %s", output_dir)
    logger.info("=" * 70)

    with ThreadPoolExecutor(max_workers=args.parallel) as executor:
        futures = {executor.submit(run_aucell, f, cfg, logger): f for f in loom_files}
        for fut in as_completed(futures):
            loom = futures[fut]
            try:
                logger.info("result: %s", fut.result())
            except Exception as exc:
                logger.error("⚠️ %s raised: %s", loom, exc)

    n_done = len([f for f in os.listdir(output_dir) if f.startswith("auc_") and f.endswith(".tsv")])
    logger.info("🎉 AUCell stage finished -- %d/%d auc_*.tsv present", n_done, len(loom_files))
    return 0


if __name__ == "__main__":
    sys.exit(main())
