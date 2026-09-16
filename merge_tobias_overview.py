"""Merge TOBIAS footprint overview.txt across three brain regions into one Parquet.

Reads <region>_Jaspar26/<celltype>_MC_footprints.bw/<motif>/<motif>_overview.txt
for AMY / HIP / PFC, keeps rows where MC_bound==1 OR MW_bound==1, and writes
a single snappy-compressed Parquet with three identifier columns appended:
brain_region, celltype, motif.

Two-phase: (1) N parallel writers each produce a chunk parquet in /tmp from a
slice of jobs; (2) chunks are concatenated into one final parquet, copied to
NFS. Phase 1 removes the single-writer bottleneck; phase 2 keeps the final
artifact single-file as requested.
"""

import glob
import os
import shutil
import time
from concurrent.futures import ThreadPoolExecutor

import pyarrow as pa
import pyarrow.compute as pacompute
import pyarrow.csv as pacsv
import pyarrow.parquet as pq


ROOTS = {
    "AMY": "/data2st2/junyi/output/atac1112/tobiasbam/AMY_Jaspar26",
    "HIP": "/data2st2/junyi/output/atac1112/tobiasbam/HIP_Jaspar26",
    "PFC": "/data2st2/junyi/output/atac1112/tobiasbam/PFC_Jaspar26",
}

OUT = "/data2st2/junyi/output/atac1112/tobiasbam/tobiasbam_all_3regions.parquet"

CHUNK_DIR = "/tmp/tobias_chunks"
N_CHUNKS = 8  # parallel parquet writers

COLS = [
    "TFBS_chr", "TFBS_start", "TFBS_end",
    "TFBS_name", "TFBS_score", "TFBS_strand",
    "peak_chr", "peak_start", "peak_end",
    "MC_score", "MW_score", "MC_bound", "MW_bound", "MC_MW_log2fc",
]

DTYPES = {
    "TFBS_chr": pa.string(),
    "TFBS_start": pa.int32(),
    "TFBS_end": pa.int32(),
    "TFBS_name": pa.string(),
    "TFBS_score": pa.float32(),
    "TFBS_strand": pa.string(),
    "peak_chr": pa.string(),
    "peak_start": pa.int32(),
    "peak_end": pa.int32(),
    "MC_score": pa.float32(),
    "MW_score": pa.float32(),
    "MC_bound": pa.int8(),
    "MW_bound": pa.int8(),
    "MC_MW_log2fc": pa.float32(),
}


def list_jobs():
    """Enumerate (brain_region, celltype, motif, txt_path) for every file."""
    jobs = []
    for br, root in ROOTS.items():
        for ct_dir_name in sorted(os.listdir(root)):
            if not ct_dir_name.endswith("_MC_footprints.bw"):
                continue
            ct_dir = os.path.join(root, ct_dir_name)
            if not os.path.isdir(ct_dir):
                continue
            celltype = ct_dir_name[: -len("_footprints.bw")]
            for motif in sorted(os.listdir(ct_dir)):
                motif_dir = os.path.join(ct_dir, motif)
                if not os.path.isdir(motif_dir):
                    continue
                txt = os.path.join(motif_dir, f"{motif}_overview.txt")
                if os.path.exists(txt):
                    jobs.append((br, celltype, motif, txt))
    return jobs


def read_one(job):
    """Read one overview.txt, drop unbound rows, append identifier columns."""
    br, celltype, motif, path = job
    tbl = pacsv.read_csv(
        path,
        read_options=pacsv.ReadOptions(column_names=COLS, skip_rows=1),
        parse_options=pacsv.ParseOptions(delimiter="\t"),
        convert_options=pacsv.ConvertOptions(column_types=DTYPES),
    )
    mask = pacompute.or_(
        pacompute.equal(tbl["MC_bound"], 1),
        pacompute.equal(tbl["MW_bound"], 1),
    )
    tbl = tbl.filter(mask)
    n = tbl.num_rows
    if n == 0:
        return None
    tbl = tbl.append_column("brain_region", pa.array([br] * n, type=pa.string()))
    tbl = tbl.append_column("celltype", pa.array([celltype] * n, type=pa.string()))
    tbl = tbl.append_column("motif", pa.array([motif] * n, type=pa.string()))
    return tbl


def write_chunk(chunk_id, job_slice):
    """One parallel writer: serialize its slice of jobs to a chunk parquet."""
    chunk_path = os.path.join(CHUNK_DIR, f"chunk_{chunk_id:02d}.parquet")
    if os.path.exists(chunk_path):
        os.remove(chunk_path)
    writer = None
    rows = 0
    written = 0
    skipped = 0
    t0 = time.time()

    for tbl in (read_one(j) for j in job_slice):
        if tbl is None:
            skipped += 1
            continue
        if writer is None:
            writer = pq.ParquetWriter(
                chunk_path,
                tbl.schema,
                compression="snappy",
                use_dictionary=True,
                write_statistics=True,
            )
        writer.write_table(tbl)
        written += 1
        rows += tbl.num_rows

    if writer is not None:
        writer.close()

    elapsed = time.time() - t0
    size_mb = os.path.getsize(chunk_path) / 1024 / 1024 if os.path.exists(chunk_path) else 0
    print(f"  chunk {chunk_id:02d}: files={written}  rows={rows:,}  "
          f"size={size_mb:.0f}MB  elapsed={elapsed:.1f}s", flush=True)
    return chunk_path


def concat_chunks(chunk_paths, out_path):
    """Concatenate per-chunk parquets into one final parquet.

    Uses parquet metadata + row-group copying to avoid re-encoding: open each
    chunk as a ParquetFile, write its row groups straight into a single
    ParquetWriter. This is O(rows) for IO, not O(rows) for re-compression.
    """
    final_tmp = "/tmp/tobiasbam_all_3regions.parquet"
    if os.path.exists(final_tmp):
        os.remove(final_tmp)

    schema = pq.read_schema(chunk_paths[0])
    writer = pq.ParquetWriter(
        final_tmp, schema,
        compression="snappy",
        use_dictionary=True,
        write_statistics=True,
    )
    t0 = time.time()
    rows = 0
    for p in chunk_paths:
        pf = pq.ParquetFile(p)
        for rg_idx in range(pf.num_row_groups):
            rg = pf.read_row_group(rg_idx)
            writer.write_table(rg)
            rows += rg.num_rows
    writer.close()
    elapsed = time.time() - t0
    print(f"Concatenated {len(chunk_paths)} chunks  rows={rows:,}  "
          f"elapsed={elapsed:.1f}s", flush=True)

    # Copy to NFS.
    print(f"Copying to {OUT} ...", flush=True)
    t1 = time.time()
    shutil.copyfile(final_tmp, out_path)
    print(f"Copied  size={os.path.getsize(out_path)/1024**3:.2f}GB  "
          f"copy_time={time.time()-t1:.1f}s", flush=True)
    os.remove(final_tmp)


def main():
    os.makedirs(CHUNK_DIR, exist_ok=True)

    t0 = time.time()
    jobs = list_jobs()
    print(f"Total files: {len(jobs)}  elapsed_list={time.time()-t0:.1f}s",
          flush=True)
    if not jobs:
        return

    # Split jobs into N chunks.
    chunk_size = (len(jobs) + N_CHUNKS - 1) // N_CHUNKS
    slices = [jobs[i * chunk_size: (i + 1) * chunk_size]
              for i in range(N_CHUNKS)]

    # Phase 1: parallel writers.
    print(f"Phase 1: {N_CHUNKS} parallel writers, {chunk_size} jobs each",
          flush=True)
    t1 = time.time()
    with ThreadPoolExecutor(max_workers=N_CHUNKS) as ex:
        chunk_paths = list(ex.map(write_chunk, range(N_CHUNKS), slices))
    print(f"Phase 1 done in {time.time()-t1:.1f}s", flush=True)

    # Phase 2: concatenate.
    print("Phase 2: concatenate chunks", flush=True)
    concat_chunks(chunk_paths, OUT)

    # Cleanup chunks.
    for p in chunk_paths:
        if os.path.exists(p):
            os.remove(p)
    try:
        os.rmdir(CHUNK_DIR)
    except OSError:
        pass

    print(f"All done in {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()