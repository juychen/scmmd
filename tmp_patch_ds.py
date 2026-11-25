#!/usr/bin/env python
"""Point the /dev/shm copy of _dataset.h5ads at /dev/shm copies of the
per-sample h5ad files (avoids NFS chunked-read latency)."""
import os
import h5py

P = "/dev/shm/atac_ds/_dataset.h5ads"
OLD = "/data2st1/junyi/output/atac0627/doublet_filtered.h5ads/anndatas"
NEW = "/dev/shm/atac_ds/anndatas"

with h5py.File(P, "r+") as f:
    g = f["uns/AnnDataSet"]
    old = [x.decode() for x in g["file_path"][...]]
    new = [x.replace(OLD, NEW) for x in old]
    g["file_path"][...] = [x.encode() for x in new]
    back = [x.decode() for x in g["file_path"][...]]

missing = [x for x in back if not os.path.exists(x)]
print("patched", len(back), "missing", len(missing))
for x in missing:
    print("MISSING", x)
for x in back[:3]:
    print("->", x)
