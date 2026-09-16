#!/usr/bin/env python3
"""Generate per-L1-celltype allc file lists from the CEMBA metadata.

Each cell in the cell-level metadata is matched to its GEO sample row via
``cell`` == ``Sample_title``; the tail of ``Sample_supplementary_file_1`` is the
allc file name.  Only files that actually exist on disk are kept, so the
resulting lists can be fed straight to ``allcools merge-allc --allc_paths``.

One list is written per ``celltype.L1`` value, e.g. ``HPF_Glut.txt``, with one
absolute path per line.
"""

from __future__ import annotations

import os
import sys

import pandas as pd

META = (
    "/data2st1/junyi/methlyatlas/mCseq/CEMBA.mC.Metadata/"
    "CEMBA.fix.mC.Metadata.selected.ictx_conflict.filtered.csv"
)
GEO = (
    "/data2st1/junyi/methlyatlas/mCseq/CEMBA.mC.Metadata/"
    "TotalGEOMetadata.2021+2023.csv"
)
ALLCS = "/data2st1/junyi/methlyatlas/mCseq/allcfiles/"
OUTDIR = "/data2st1/junyi/methlyatlas/mCseq/CEMBA.mC.Metadata/"


def main() -> int:
    meta = pd.read_csv(META, usecols=["cell", "celltype.L1"])
    geo = pd.read_csv(
        GEO,
        usecols=["Sample_title", "Sample_supplementary_file_1"],
        low_memory=False,
    )

    df = meta.merge(geo, left_on="cell", right_on="Sample_title", how="left")
    n_unmapped = int(df["Sample_supplementary_file_1"].isna().sum())

    df["allc_path"] = ALLCS + df["Sample_supplementary_file_1"].str.split("/").str[-1]

    def _exists(path):
        return isinstance(path, str) and os.path.exists(path)

    exists = df["allc_path"].map(_exists)
    n_missing = int((df["allc_path"].notna() & ~exists).sum())
    df = df[exists & df["celltype.L1"].notna()]

    summary = []
    for l1, sub in df.groupby("celltype.L1"):
        paths = sorted(set(sub["allc_path"]))
        if not paths:
            continue
        with open(os.path.join(OUTDIR, f"{l1}.txt"), "w") as f:
            f.write("\n".join(paths) + "\n")
        summary.append((l1, len(paths)))

    total = sum(n for _, n in summary)
    print(f"cells in metadata: {len(meta):,}")
    print(f"cells without a GEO row: {n_unmapped:,}")
    print(f"cells whose allc file is absent on disk: {n_missing:,}")
    print(f"wrote {len(summary)} lists, {total:,} paths total -> {OUTDIR}")
    for l1, n in sorted(summary, key=lambda x: (-x[1], x[0])):
        print(f"  {l1:20s} {n:,}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
