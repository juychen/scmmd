#!/usr/bin/env python3
"""Build per-(region x condition x celltype.L1_ct) ATAC bigWigs, mmllm-style.

Reuses the existing TOBIAS per-sample per-L2-cluster BAMs:
    /data2st2/junyi/output/atac1112/tobiasbam/<sample>/<region>_<L2>_<cond>.bam
merges them per L1_ct, then runs the exact bamCoverage command used for mmllm's
6 bulk tracks (CPM, binSize 10, mm10 blacklist, extend/center reads, dedup, MAPQ30).

Merged BAMs are ~= the size of their inputs (Glut in AMY_MC is ~106 GB), so they
are written to the first temp location that has room:
    1. /dev/shm/l1_tracks_tmp        (RAM disk, 756 GB, fastest)
    2. /data1st2/junyi/tmp_l1tracks  (disk, 5.7 TB free)
Never to /tmp or to /data2st1 (both nearly full).

Usage:
    python build_l1_ct_tracks.py --all --jobs 3
    python build_l1_ct_tracks.py --regions AMY --conditions MC --l1cts Glut
    python build_l1_ct_tracks.py --all --dry-run
"""
import argparse
import os
import shutil
import subprocess
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = "/data2st2/junyi/output/atac1112/tobiasbam"
OUT = "/data2st1/junyi/output/atac1112/l1_tracks"
TMP_CANDIDATES = [
    "/dev/shm/l1_tracks_tmp",
    "/data1st2/junyi/tmp_l1tracks",
]
BLACKLIST = "/data2st1/junyi/ref/mm10-blacklist.v2.bed"
BAMCOV = "/home/junyichen/anaconda3/envs/tobias/bin/bamCoverage"
M5C = "/data2st1/junyi/output/atac1112/3REGIONS_5mc_new.h5ads"

REGIONS = ["AMY", "HIP", "PFC"]
CONDS = ["MC", "MW"]
L1CTS = ["Glut", "GABA", "Oligo", "Microglia", "Astrocyte", "OPC", "Vascular", "Epen", "Immune"]

GI = 1024 ** 3
_lock = threading.Lock()
_reserved = {p: 0 for p in TMP_CANDIDATES}


def log(msg):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def pick_tmp(need_bytes):
    """Reserve space in the first temp location that can hold this merge."""
    with _lock:
        for p in TMP_CANDIDATES:
            os.makedirs(p, exist_ok=True)
            free = shutil.disk_usage(p).free
            if free - _reserved[p] > need_bytes + 2 * GI:
                _reserved[p] += need_bytes
                return p
    return None


def release_tmp(p, need_bytes):
    with _lock:
        _reserved[p] = max(0, _reserved[p] - need_bytes)


def run(cmd, tag):
    log(f"RUN {tag}")
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        log(f"FAIL {tag} rc={p.returncode}")
        log(((p.stdout or "") + (p.stderr or ""))[-2000:])
        return False
    log(f"OK {tag} in {time.time() - t0:.0f}s")
    return True


def collect_inputs(obs, region, cond, l1ct):
    sub = obs[(obs["Region"] == region) & (obs["Condition"] == cond) & (obs["celltype.L1_ct"] == l1ct)]
    if len(sub) == 0:
        return None, 0, []
    files = []
    for s in sorted(sub["Sample_name"].unique()):
        for l2 in sorted(sub.loc[sub["Sample_name"] == s, "celltype.L2"].unique()):
            f = os.path.join(ROOT, s, f"{region}_{l2}_{cond}.bam")
            if os.path.exists(f):
                files.append(f)
    return sub, len(sub), files


def build_one(obs, region, cond, l1ct):
    tag = f"{region}_{cond}_{l1ct}"
    bw = os.path.join(OUT, f"{tag}_track.bw")
    if os.path.exists(bw) and os.path.getsize(bw) > 10 * 1024 * 1024:
        log(f"SKIP {tag}: already exists ({os.path.getsize(bw)/1e6:.0f} MB)")
        return "skip"

    sub, ncells, files = collect_inputs(obs, region, cond, l1ct)
    if not files:
        log(f"SKIP {tag}: no input bams (cells={ncells})")
        return "skip"

    need = int(sum(os.path.getsize(f) for f in files) * 1.05)
    log(f"== {tag}: {ncells} cells, {len(files)} bams, {need/GI:.1f} GB")

    tmpdir = pick_tmp(need)
    if tmpdir is None:
        log(f"DEFER {tag}: no temp space for {need/GI:.1f} GB")
        return "defer"

    merged = os.path.join(tmpdir, f"{tag}.bam")
    try:
        if not run(["samtools", "merge", "-f", "-@", "16", merged] + files, f"merge {tag}"):
            return "fail"
        if not run(["samtools", "index", "-@", "8", merged], f"index {tag}"):
            return "fail"
        cmd = [
            BAMCOV, "-b", merged, "-o", bw,
            "--blackListFileName", BLACKLIST,
            "--binSize", "10",
            "--normalizeUsing", "CPM",
            "--effectiveGenomeSize", "2652783500",
            "--numberOfProcessors", "32",
            "--extendReads", "--centerReads", "--ignoreDuplicates",
            "--minMappingQuality", "30",
        ]
        if not run(cmd, f"bamCoverage {tag}"):
            if os.path.exists(bw):
                os.remove(bw)
            return "fail"
        log(f"DONE {tag}: {os.path.getsize(bw)/1e6:.0f} MB")
        return "ok"
    finally:
        for f in (merged, merged + ".bai"):
            if os.path.exists(f):
                try:
                    os.remove(f)
                except OSError:
                    pass
        release_tmp(tmpdir, need)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--regions", nargs="+", default=REGIONS)
    ap.add_argument("--conditions", nargs="+", default=CONDS)
    ap.add_argument("--l1cts", nargs="+", default=L1CTS)
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    log(f"blacklist={os.path.exists(BLACKLIST)} bamCoverage={os.path.exists(BAMCOV)}")
    for p in TMP_CANDIDATES:
        os.makedirs(p, exist_ok=True)
        log(f"tmp {p}: {shutil.disk_usage(p).free/GI:.0f} GB free")

    import anndata as ad

    a = ad.read_h5ad(M5C, backed="r")
    obs = a.obs[["Sample_name", "Region", "Condition", "celltype.L2", "celltype.L1_ct"]].copy()
    log(f"obs cells: {len(obs)}")

    combos = [(r, c, k) for r in args.regions for c in args.conditions for k in args.l1cts]
    log(f"planned {len(combos)} tracks, jobs={args.jobs}")

    if args.dry_run:
        for r, c, k in combos:
            _, n, files = collect_inputs(obs, r, c, k)
            gb = sum(os.path.getsize(f) for f in files) / GI if files else 0
            exists = os.path.exists(os.path.join(OUT, f"{r}_{c}_{k}_track.bw"))
            log(f"  {r}_{c}_{k}: cells={n} bams={len(files)} {gb:.1f} GB "
                f"{'EXISTS' if exists else ''}")
        return

    stats = {"ok": 0, "skip": 0, "fail": 0, "defer": 0}
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(build_one, obs, r, c, k): f"{r}_{c}_{k}" for r, c, k in combos}
        for fut in as_completed(futs):
            tag = futs[fut]
            try:
                stats[fut.result()] += 1
            except Exception as e:
                log(f"EXC {tag}: {e!r}")
                stats["fail"] += 1
    log(f"ALL DONE {stats}")


if __name__ == "__main__":
    main()
