#!/usr/bin/env python
"""Compare snapatac2 RPKM/bin10 tracks against the bamCoverage CPM tracks.

usage: verify_tracks.py <snap.bw> <ref.bw> [chrom ...]
"""
import sys

import numpy as np
import pyBigWig


def load(bw, chrom, step=1):
    n = bw.chroms()[chrom]
    v = np.nan_to_num(bw.values(chrom, 0, n, numpy=True), nan=0.0)
    return v.astype(np.float64)


def main():
    a, b = sys.argv[1], sys.argv[2]
    chroms = sys.argv[3:] or ["chr1", "chr2"]
    ba, bb = pyBigWig.open(a), pyBigWig.open(b)
    print(f"snap={a}\nref ={b}")
    print("chroms snap:", len(ba.chroms()), "ref:", len(bb.chroms()))
    for c in chroms:
        x = load(ba, c)
        y = load(bb, c)
        if x.shape != y.shape:
            m = min(x.shape[0], y.shape[0])
            x, y = x[:m], y[:m]
        nz = (x > 0) | (y > 0)
        r = np.corrcoef(x[nz], y[nz])[0, 1]
        sp = np.corrcoef(np.argsort(np.argsort(x[nz])),
                         np.argsort(np.argsort(y[nz])))[0, 1]
        print(f"{c}: n={nz.sum()} pearson={r:.4f} spearman~{sp:.4f} "
              f"sum_snap={x.sum():.3g} sum_ref={y.sum():.3g} "
              f"ratio={x.sum() / y.sum():.3f} mean_snap={x.mean():.3g} "
              f"mean_ref={y.mean():.3g} ratio_mean={x.mean() / y.mean():.3f}")
    for pos in [(10000000, 10000010), (30000000, 30000010)]:
        print("chr1", pos, "snap", ba.values("chr1", pos[0], pos[1]),
              "ref", bb.values("chr1", pos[0], pos[1]))


if __name__ == "__main__":
    main()
