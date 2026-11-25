#!/usr/bin/env python
"""Export CPM bigWig tracks with SnapATAC2, replicating the bamCoverage
settings used for the BULK group tracks (see 17trackPPbulk.sh):

    bamCoverage --binSize 10 --normalizeUsing CPM \
        --blackListFileName /data2st1/junyi/ref/mm10-blacklist.v2.bed \
        --extendReads --centerReads --ignoreDuplicates --minMappingQuality 30

SnapATAC2 equivalent:
    snap.ex.export_coverage(..., bin_size=10, normalization="RPKM",
        counting_strategy="fragment", blacklist=mm10-blacklist.v2.bed,
        max_frag_length=2000)

Note on scale: snapatac2's "RPKM" divides the per-bin count by
(sum_fragment_bases * bin_size / 1e9); with bin_size=10 that equals
count/(sum_bases/1e8) = 100 * its "CPM", which is numerically on the
same scale as bamCoverage --normalizeUsing CPM (verified r~0.98,
slope ~1.05 on chr1/chr2 vs the bamCoverage AMY_MC_Glut pilot).

Levels
------
sample     : REGION_COND  (6 tracks, e.g. AMY_MC_track.bw)
sample_l1  : REGION_COND_L1ct (54 tracks, e.g. AMY_MC_Glut_track.bw)
l1         : L1ct (9 tracks, e.g. Glut_track.bw)

Usage
-----
python 19export_bw_snapatac.py --levels sample sample_l1 l1
python 19export_bw_snapatac.py --levels sample_l1 --only AMY_MC HIP_MW
"""
import argparse
import os
import sys
import time

import numpy as np
import snapatac2 as snap

DS = "/data2st1/junyi/output/atac0627/doublet_filtered.h5ads/_dataset.h5ads"
BLACKLIST = "/data2st1/junyi/ref/mm10-blacklist.v2.bed"
OUT_ROOT = "/data2st1/junyi/output/atac1112/tracks_snapatac"
TMPDIR = "/dev/shm/snapatac_tmp"

L1_ORDER = ["Glut", "GABA", "Oligo", "Microglia", "Astrocyte",
            "OPC", "Vascular", "Epen", "Immune"]


def build_labels(ds, level):
    reg = np.asarray(ds.obs["Region"])
    con = np.asarray(ds.obs["Condition"])
    l1 = np.asarray(ds.obs["celltype.L1_ct"])
    if level == "sample":
        return np.array([f"{r}_{c}" for r, c in zip(reg, con)], dtype=object)
    if level == "sample_l1":
        return np.array([f"{r}_{c}_{t}" for r, c, t in zip(reg, con, l1)],
                        dtype=object)
    if level == "l1":
        return np.asarray(l1, dtype=object)
    raise ValueError(level)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", nargs="+", default=["sample", "sample_l1", "l1"],
                    choices=["sample", "sample_l1", "l1"])
    ap.add_argument("--only", nargs="*", default=None,
                    help="only export these group names (optional)")
    ap.add_argument("--dataset", default=DS,
                    help="AnnDataSet path (default: original on /data2st1)")
    ap.add_argument("--out-root", default=OUT_ROOT)
    ap.add_argument("--jobs", type=int, default=32)
    args = ap.parse_args()

    t0 = time.time()
    print(f"[{time.strftime('%F %T')}] reading {args.dataset}", flush=True)
    ds = snap.read_dataset(args.dataset)
    print(f"loaded {ds.shape} in {time.time() - t0:.1f}s", flush=True)

    for level in args.levels:
        labels = build_labels(ds, level)
        groups = sorted(set(labels.tolist()))
        if args.only:
            groups = [g for g in groups
                      if any(g == o or g.startswith(o + "_") for o in args.only)]
        if not groups:
            print(f"level {level}: no groups selected, skip", flush=True)
            continue
        outdir = os.path.join(args.out_root, level)
        os.makedirs(outdir, exist_ok=True)
        os.makedirs(TMPDIR, exist_ok=True)
        for g in groups:
            n = int(np.sum(labels == g))
            print(f"[{time.strftime('%F %T')}] export {level}/{g}: {n} cells",
                  flush=True)
        t = time.time()
        res = snap.ex.export_coverage(
            ds, groupby=labels.tolist(), selections=groups,
            bin_size=10, normalization="RPKM", blacklist=BLACKLIST,
            max_frag_length=2000, counting_strategy="fragment",
            out_dir=outdir, prefix="", suffix="_track.bw",
            tempdir=TMPDIR, n_jobs=args.jobs,
        )
        for k, v in res.items():
            print(f"  -> {k}: {v} ({os.path.getsize(v) / 1e9:.3f} GB)", flush=True)
        print(f"[{time.strftime('%F %T')}] level {level} done in "
              f"{time.time() - t:.1f}s", flush=True)

    print(f"ALL DONE in {time.time() - t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
