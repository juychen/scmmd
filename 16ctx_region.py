#!/usr/bin/env python3
"""PySCENIC CTX (cisTarget regulon pruning) for the region-level looms.

Same scheme as `ctx_sn0827.py`:
    pyscenic ctx adj_<region>.tsv <rankings.feather>
        --annotations_fname <motifs.tbl>
        --expression_mtx_fname <region>.loom
        --output ctx_<region>.csv --num_workers N

Defaults are the v9nr resources; pass --annotations/--database (and
--output-dir) to run against another database, e.g. the mouse v10nr_clust set:

    python 16ctx_region.py \
        --annotations /data2st1/junyi/scenic/mouse/motif/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl \
        --database    /data2st1/junyi/scenic/mouse/motif/mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
        --output-dir  /data1st1/junyi/output/sn1029Region/scenic_ctx_v10

Output : <output-dir>/ctx_<region>.csv   (default <input-dir>/scenic_ctx)
Logs   : <output-dir>/ctx_run.log and logs/<region>.ctx.log
"""

import argparse
import logging
import multiprocessing
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd

# === DEFAULT PATHS (v9nr) ===
INPUT_DIR = "/data1st1/junyi/output/sn1029Region"
ANNOTATIONS = "/data1st2/yejun/scenic/mouse/motifs-v9-nr.mgi-m0.001-o0.0.tbl"
DATABASE = ("/data1st2/yejun/scenic/mouse/"
            "mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather")

DEFAULT_PARALLEL = 4


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="PySCENIC ctx on region-level looms")
    p.add_argument("--input-dir", default=INPUT_DIR, help="directory with the .loom files")
    p.add_argument("--grn-dir", default=None, help="directory with adj_*.tsv (default <input-dir>/scenic_grn)")
    p.add_argument("--output-dir", default=None, help="ctx output directory (default <input-dir>/scenic_ctx)")
    p.add_argument("--annotations", default=ANNOTATIONS, help="motif annotation .tbl")
    p.add_argument("--database", action="append", default=None,
                   help="cisTarget ranking .feather, repeatable (default: the v9nr 10kb database)")
    p.add_argument("--parallel", type=int, default=DEFAULT_PARALLEL, help="looms processed at once")
    p.add_argument("--workers", type=int, default=None, help="--num_workers (default: ncpu-4, max 32)")
    p.add_argument("--regions", default=None, help="comma separated subset, e.g. AMY,STR")
    return p.parse_args(argv)


def fix_grn_columns(grn_path, logger):
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


def run_ctx(loom_file, cfg, logger):
    base = os.path.splitext(loom_file)[0]
    expr_path = os.path.join(cfg["input_dir"], loom_file)
    grn_path = os.path.join(cfg["grn_dir"], f"adj_{base}.tsv")
    output_path = os.path.join(cfg["output_dir"], f"ctx_{base}.csv")
    flag_file = os.path.join(cfg["output_dir"], f"{base}.processing")
    job_log = os.path.join(cfg["log_dir"], f"{base}.ctx.log")

    if os.path.exists(flag_file):
        logger.info("⏩ skip (processing flag): %s", base)
        return f"skip-flag: {base}"
    if os.path.exists(output_path) and os.path.getsize(output_path) > 1024:
        logger.info("✅ skip (output exists): %s", base)
        return f"skip-exists: {base}"
    if not os.path.exists(grn_path) or os.path.getsize(grn_path) < 1024:
        logger.warning("🚫 skip (GRN missing/empty): %s", base)
        return f"missing-grn: {base}"

    with open(flag_file, "w") as fh:
        fh.write(f"Processing {loom_file}\n")

    if not fix_grn_columns(grn_path, logger):
        os.remove(flag_file)
        return f"bad-grn: {base}"

    cmd = ["pyscenic", "ctx", grn_path, *cfg["databases"],
           "--annotations_fname", cfg["annotations"],
           "--expression_mtx_fname", expr_path,
           "--output", output_path,
           "--num_workers", str(cfg["workers"])]

    t0 = time.time()
    logger.info("▶ CTX start: %s (workers=%d, db=%s)", base, cfg["workers"],
                os.path.basename(cfg["databases"][0]))
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


def main(argv=None):
    args = parse_args(argv)

    input_dir = args.input_dir
    grn_dir = args.grn_dir or os.path.join(input_dir, "scenic_grn")
    output_dir = args.output_dir or os.path.join(input_dir, "scenic_ctx")
    log_dir = os.path.join(output_dir, "logs")
    databases = args.database or [DATABASE]
    workers = args.workers or min(32, max(1, multiprocessing.cpu_count() - 4))

    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stdout),
                  logging.FileHandler(os.path.join(output_dir, "ctx_run.log"))],
    )
    logger = logging.getLogger(__name__)

    loom_files = sorted(f for f in os.listdir(input_dir) if f.endswith(".loom"))
    if args.regions:
        wanted = {r.strip() for r in args.regions.split(",") if r.strip()}
        loom_files = [f for f in loom_files if os.path.splitext(f)[0] in wanted]
    if not loom_files:
        logger.error("no .loom in %s", input_dir)
        return 1

    cfg = dict(input_dir=input_dir, grn_dir=grn_dir, output_dir=output_dir,
               log_dir=log_dir, annotations=args.annotations,
               databases=databases, workers=workers)

    logger.info("=" * 70)
    logger.info("PySCENIC CTX -- %d looms", len(loom_files))
    logger.info("parallel jobs=%d, --num_workers=%d", args.parallel, workers)
    logger.info("annotations: %s", args.annotations)
    for db in databases:
        logger.info("database   : %s", db)
    logger.info("output     : %s", output_dir)
    logger.info("=" * 70)

    with ThreadPoolExecutor(max_workers=args.parallel) as executor:
        futures = {executor.submit(run_ctx, f, cfg, logger): f for f in loom_files}
        for fut in as_completed(futures):
            loom = futures[fut]
            try:
                logger.info("result: %s", fut.result())
            except Exception as exc:
                logger.error("⚠️ %s raised: %s", loom, exc)

    n_done = len([f for f in os.listdir(output_dir) if f.startswith("ctx_") and f.endswith(".csv")])
    logger.info("🎉 CTX stage finished -- %d/%d ctx_*.csv present", n_done, len(loom_files))
    return 0


if __name__ == "__main__":
    sys.exit(main())
