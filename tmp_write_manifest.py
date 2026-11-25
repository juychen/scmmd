#!/usr/bin/env python
"""Write README + manifest for the snapatac2 bigWig tracks."""
import glob
import json
import os

import numpy as np
import snapatac2 as snap

OUT = "/data2st1/junyi/output/atac1112/tracks_snapatac"
DS = "/dev/shm/atac_ds/_dataset.h5ads"
TINY = 50

ds = snap.read_dataset(DS)
reg = np.asarray(ds.obs["Region"])
con = np.asarray(ds.obs["Condition"])
l1 = np.asarray(ds.obs["celltype.L1_ct"])
samp = np.array([f"{r}_{c}" for r, c in zip(reg, con)])
lab = np.array([f"{r}_{c}_{t}" for r, c, t in zip(reg, con, l1)])

rows = []
for level, labv in [("sample", samp), ("sample_l1", lab)]:
    for p in sorted(glob.glob(f"{OUT}/{level}/*_track.bw")):
        g = os.path.basename(p)[:-len("_track.bw")]
        n = int((labv == g).sum())
        rows.append((level, g, n, os.path.getsize(p) / 1e9, p))

with open(f"{OUT}/manifest.tsv", "w") as fh:
    fh.write("level\tgroup\tcells\tsize_GB\tpath\n")
    for r in rows:
        fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.3f}\t{r[4]}\n")

tiny = [r for r in rows if r[2] < TINY]
small = [r for r in rows if r[2] < 200]

readme = f"""# ATAC bigWig tracks produced with SnapATAC2 (mm10, 10-bp bins)

Generated {__import__('time').strftime('%F %T')} from
`/data2st1/junyi/output/atac0627/doublet_filtered.h5ads/_dataset.h5ads`
(176,318 QC-passed cells; fragments copied to `/dev/shm/atac_ds` for speed).

## Why RPKM instead of CPM

`bamCoverage --normalizeUsing CPM --binSize 10` (the mmllm BULK setting in
`17trackPPbulk.sh:37-48`) is numerically reproduced by SnapATAC2
`normalization="RPKM"` with `bin_size=10`. In SnapATAC2's implementation
(`snapatac2-core/src/export.rs`) the norm factor is the **total fragment
bases**, and `RPKM = count / (bases * bin_size / 1e9)`, so
`RPKM = CPM_snap * bin_size / 1000 = CPM_snap * 100` for 10-bp bins.
Using `normalization="CPM"` would give values 100x smaller than bamCoverage.

## Settings (copied from 17trackPPbulk.sh:37-48)

| bamCoverage | SnapATAC2 equivalent used here |
|---|---|
| `--binSize 10` | `bin_size=10` |
| `--normalizeUsing CPM` | `normalization="RPKM"` (see above) |
| `--blackListFileName mm10-blacklist.v2.bed` | `blacklist="/data2st1/junyi/ref/mm10-blacklist.v2.bed"` |
| fragment counting (`--extendReads --centerReads`) | `counting_strategy="fragment"` |
| fragment length cap | `max_frag_length=2000` |
| `--numberOfProcessors 32` | `n_jobs=32` |

Not reproducible from the fragment matrix: `--ignoreDuplicates` and
`--minMappingQuality 30` (SnapATAC2 stores already-filtered fragments; no
per-read MAPQ/dedup filters are applied at export time).

## Outputs

* `sample/`  – 6 tracks `{{REGION}}_{{COND}}_track.bw`
  (AMY/HIP/PFC x MC/MW). Cell numbers: AMY_MC 25,400; AMY_MW 27,381;
  HIP_MC 23,354; HIP_MW 33,589; PFC_MC 32,360; PFC_MW 34,234.
* `sample_l1/` – {sum(1 for r in rows if r[0] == 'sample_l1')} tracks
  `{{REGION}}_{{COND}}_{{L1ct}}_track.bw` with `celltype.L1_ct` in
  Glut, GABA, Oligo, Microglia, Astrocyte, OPC, Vascular, Epen, Immune.
  Only combinations that actually contain cells were exported
  (HIP_MW_Immune has 0 cells).
* `manifest.tsv` – group, cell number, size, path.

## Validation vs bamCoverage

| comparison | Pearson (chr1) | scale snap/ref |
|---|---|---|
| AMY_MC sample vs BULK `AMY_MC_track.bw` | 0.968 | 0.93 |
| AMY_MC_Glut vs bamCoverage L1 pilot `l1_tracks/AMY_MC_Glut_track.bw` | 0.983 | 0.95 |
| AMY_MW_Glut vs BULK `AMY_MW_track.bw` | 0.965 | 1.20 (Glut vs all cells) |

Shape agreement is ~0.96-0.98; the residual ~5% scale offset comes from
SnapATAC2 dividing by total fragment bases instead of read count, and from
the slightly different cell sets (BULK includes pre-QC cells).

## Caveat: very small groups

These tracks are normalized **per track**, so groups with very few cells get
inflated values (e.g. PFC_MC_Immune, 2 cells, values ~50-160 instead of ~0.03):

{chr(10).join(f'* {r[1]}: {r[2]} cells ({r[3]*1e3:.1f} MB)' for r in tiny)}

Groups with <200 cells: {', '.join(f'{r[1]} ({r[2]})' for r in small)}.
Interpret small-group tracks qualitatively only (or pool them).

## Reproduce

```bash
python /home/junyichen/code/scmmd/19export_bw_snapatac.py \\
    --levels sample sample_l1 \\
    --dataset /data2st1/junyi/output/atac0627/doublet_filtered.h5ads/_dataset.h5ads \\
    --out-root /data2st1/junyi/output/atac1112/tracks_snapatac --jobs 32
```

The `/dev/shm/atac_ds` copy (43 GB, tmpfs) only exists in RAM; rebuild it with
`cp -a` and `tmp_patch_ds.py` (see /home/junyichen/code/scmmd), or just run
against the original dataset path (slower: ~6-12 MB/s NFS chunk reads).
Remove it with `rm -rf /dev/shm/atac_ds /dev/shm/snapatac_tmp` when done.
"""

with open(f"{OUT}/README.md", "w") as fh:
    fh.write(readme)
print("wrote", f"{OUT}/README.md", "and manifest.tsv with", len(rows), "rows")
