#!/usr/bin/env python3
"""Build ONE loom file per brain region from the CUMS / CSDS single-nucleus h5ads.

This is the "region-level" companion of ``11subsetloom_new.py`` /
``11subsetloom_3v3.py``: instead of one loom per *cell type*, it writes one loom
per *brain region*, where the cells of a loom are the union of a few selected
experimental groups (condition x sex x region).

Groups built by default (edit ``GROUP_SPECS`` below, or use ``--config`` /
``--groups`` to change them without touching the code):

    CUMS_M   CUMS  SUS   male    AMY PFC TH HIP
    CUMS_F   CUMS  SUS   female  AMY PFC TH HIP
    CSDS_M   CSDS  CSDS  male    AMY PFC TH HIP STR HY MB
    CON_M    CUMS  CON   male    AMY PFC TH HIP STR HY MB
    CON_F    CUMS  CON   female  AMY PFC TH HIP

which gives one loom per region:

    AMY / PFC / TH / HIP   <- CUMS_M + CUMS_F + CSDS_M + CON_M + CON_F
    STR / HY  / MB         <- CSDS_M + CON_M

Every loom carries the per-cell metadata needed downstream (AUCell, differential
regulon analysis): ``Group``, ``Condition``, ``Sex``, ``Dataset``, ``Region``,
``RegionSource``, ``Sample``, ``Donor``, ``celltype.L1``, ``celltype.L2``,
``nGene``, ``nUMI``.

Notes
-----
* The male control cells are physically the same cells in the CUMS and the CSDS
  h5ad (verified: all 289,670 CSDS ``CON`` barcodes are contained in the CUMS
  ``CON`` cells).  Cells are therefore de-duplicated by barcode across datasets,
  so a shared control cell is never counted twice.
* The matrix is read directly with ``h5py`` (CSR row slicing) instead of
  ``anndata.read_h5ad(backed="r")``: the obs table of the CUMS file is 1.1 M
  cells x 69 columns and the backed reader needs minutes and several GB just to
  load obs.  ``--reader anndata`` keeps a backed-reader fallback for files whose
  matrix is not stored as CSR.
* ``--layer counts`` (default) writes the raw UMI counts, i.e. the same kind of
  matrix as the existing ``sn0827L1/*.loom`` files.  Use ``--layer X`` for the
  normalised (log1p CPM) matrix.
* The loom itself is written by a small h5py writer (``create_loom_fast``)
  rather than ``loompy.create``: loompy's sparse path densifies column windows
  into (64, 64) gzip chunks, which needs millions of tiny HDF5 writes and takes
  hours for a region-scale matrix.  The output is a plain Loom 3.0 file that
  ``loompy.connect``, pySCENIC and SCopeLoomR read unchanged.
* ``--max_cells_per_group`` randomly caps the cells of each group inside each
  region (seeded), which is the knob to pull when the looms get too big for
  pySCENIC.

Examples
--------
    # only inspect the selection (no file is written)
    python 11subsetloom_meth.py --dry-run

    # build every region
    python 11subsetloom_meth.py --output /data1st1/junyi/output/meth_region_loom

    # build two regions only, 20k cells per group per region, overwrite
    python 11subsetloom_meth.py --regions AMY,STR --max_cells_per_group 20000 \
        --output /data1st1/junyi/output/meth_region_loom --overwrite
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

# --------------------------------------------------------------------------- #
# Default configuration -- edit here, or override with --config / --input
# --------------------------------------------------------------------------- #

# dataset name -> h5ad path
DEFAULT_INPUTS = {
    "CUMS": "/data2st2/junyi/output/stg1028/CUMS_4VN/CUMS_4VN.h5ad",
    "CSDS": "/data2st2/junyi/output/stg1028/CSDS_4VN/CSDS_4VN.h5ad",
}

DEFAULT_OUTPUT = "/data1st1/junyi/output/meth_region_loom"

# One entry per experimental group.  ``regions`` are values of the obs ``region``
# column (the data uses ``HPF`` for hippocampus; it is renamed on output).
GROUP_SPECS = [
    dict(dataset="CUMS", status="SUS",  label="CUMS", sex="M",
         regions=("AMY", "PFC", "TH", "HPF")),
    dict(dataset="CUMS", status="SUS",  label="CUMS", sex="F",
         regions=("AMY", "PFC", "TH", "HPF")),
    dict(dataset="CSDS", status="CSDS", label="CSDS", sex="M",
         regions=("AMY", "PFC", "TH", "HPF", "STR", "HY", "MB")),
    dict(dataset="CUMS", status="CON",  label="CON",  sex="M",
         regions=("AMY", "PFC", "TH", "HPF", "STR", "HY", "MB")),
    dict(dataset="CUMS", status="CON",  label="CON",  sex="F",
         regions=("AMY", "PFC", "TH", "HPF")),
]

# order in which the looms are written / printed
REGION_ORDER = ["AMY", "PFC", "TH", "HPF", "STR", "HY", "MB"]

# obs columns holding the metadata
STATUS_COL = "status"
SEX_COL = "sex"
REGION_COL = "region"
SAMPLE_COL = "sample"
DONOR_COL = "donor"
CELLTYPE_COLS = ("celltype.L1", "celltype.L2")

# output name of a region (HPF = hippocampal formation = "HIP" in the ATAC data)
DEFAULT_RENAME = {"HPF": "HIP"}


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

def log(msg: str = "") -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def decode_str_array(arr) -> np.ndarray:
    """bytes / vlen-str / object array -> numpy object array of python str."""
    arr = np.asarray(arr)
    if arr.dtype.kind == "S":
        return np.char.decode(arr, "utf-8")
    if arr.dtype.kind == "O":
        return np.array(
            [x.decode("utf-8") if isinstance(x, (bytes, np.bytes_)) else str(x) for x in arr],
            dtype=object,
        )
    return arr.astype(str)


class ObsColumn:
    """One ``obs`` column, kept encoded until a subset is requested.

    Categorical columns stay as ``codes`` (1-2 bytes/cell) and are only decoded
    for the selected cells, so touching the obs of a 1.1 M-cell file is cheap.
    """

    def __init__(self, group: h5py.Group | h5py.Dataset):
        self.is_cat = isinstance(group, h5py.Group)
        if self.is_cat:
            self.cats = decode_str_array(group["categories"][:])
            self.codes = group["codes"][:]
        else:
            self.values = decode_str_array(group[:])

    def match(self, values) -> np.ndarray:
        """Boolean mask of the rows whose value is in ``values``."""
        wanted = {str(v) for v in values}
        if self.is_cat:
            hit = np.array(
                [i for i, c in enumerate(self.cats) if c in wanted], dtype=self.codes.dtype
            )
            return np.isin(self.codes, hit)
        return np.isin(self.values, list(wanted))

    def take(self, idx: np.ndarray) -> np.ndarray:
        """Decoded values at the given row indices."""
        idx = np.asarray(idx, dtype=np.int64)
        if self.is_cat:
            codes = self.codes[idx]
            out = np.empty(len(idx), dtype=object)
            ok = codes >= 0
            out[ok] = self.cats[codes[ok]]
            out[~ok] = ""
            return out
        return self.values[idx]


class CsrSource:
    """Row-sliceable CSR matrix read straight from an h5ad (``X`` or a layer).

    ``data`` / ``indices`` are loaded once per file (numpy row slicing after that
    is essentially free) instead of h5py fancy-indexing every block.
    """

    def __init__(self, path: str, layer: str | None = None):
        self.path = path
        self.layer = layer
        self._h5 = h5py.File(path, "r")
        grp = self._group()
        enc = grp.attrs.get("encoding-type", "csr_matrix")
        if isinstance(enc, bytes):
            enc = enc.decode()
        if enc not in ("csr_matrix", None):
            raise ValueError(
                f"{path}: layer {layer or 'X'!r} is stored as {enc!r}; "
                "only csr_matrix is supported by the fast reader, "
                "re-run with --reader anndata"
            )
        self.shape = tuple(int(x) for x in grp.attrs["shape"])
        self._indptr = grp["indptr"][:]
        self._data = None
        self._indices = None

    def _group(self) -> h5py.Group:
        return self._h5["X"] if self.layer is None else self._h5["layers"][self.layer]

    def _ensure_loaded(self) -> None:
        if self._data is None:
            grp = self._group()
            self._data = grp["data"][:]
            self._indices = grp["indices"][:]
            log(f"    cached {os.path.basename(self.path)}[{self.layer or 'X'}]: "
                f"{len(self._data):,} non-zeros (~{len(self._data) * 8 / 1e9:.1f} GB)")

    def rows(self, idx: np.ndarray) -> sp.csr_matrix:
        """CSR sub-matrix for (ascending) row indices ``idx``."""
        self._ensure_loaded()
        idx = np.asarray(idx, dtype=np.int64)
        starts = self._indptr[idx]
        lens = (self._indptr[idx + 1] - starts).astype(np.int64)
        total = int(lens.sum())

        new_indptr = np.zeros(len(idx) + 1, dtype=np.int64)
        np.cumsum(lens, out=new_indptr[1:])
        if total:
            pos = np.repeat(starts - new_indptr[:-1], lens) + np.arange(total, dtype=np.int64)
            data = self._data[pos]
            indices = self._indices[pos]
        else:
            data = np.zeros(0, dtype=self._data.dtype)
            indices = np.zeros(0, dtype=self._indices.dtype)
        return sp.csr_matrix((data, indices, new_indptr), shape=(len(idx), self.shape[1]))

    def close(self) -> None:
        self._data = self._indices = None
        self._h5.close()


class AnndataSource:
    """Backed-anndata fallback reader (slow on obs, but works for any encoding)."""

    def __init__(self, path: str, layer: str | None = None):
        import anndata

        self.layer = layer
        self._ad = anndata.read_h5ad(path, backed="r")
        self.shape = tuple(self._ad.shape)

    def rows(self, idx: np.ndarray) -> sp.csr_matrix:
        sub = self._ad[np.asarray(idx, dtype=np.int64), :].to_memory()
        mat = sub.X if self.layer is None else sub.layers[self.layer]
        return sp.csr_matrix(mat)

    def close(self) -> None:
        try:
            self._ad.file.close()
        except Exception:
            pass


class Dataset:
    """One input h5ad: obs columns + matrix source."""

    def __init__(self, name: str, path: str, layer: str | None, reader: str = "h5py"):
        self.name = name
        self.path = path
        if not os.path.exists(path):
            raise FileNotFoundError(f"{name}: {path} does not exist")
        self._h5 = h5py.File(path, "r")
        if layer is not None and layer not in self._h5["layers"]:
            raise KeyError(
                f"{name}: no layer {layer!r}; available layers: "
                f"{sorted(self._h5['layers'].keys())}"
            )
        self.index = decode_str_array(self._h5["obs"]["_index"][:])
        self.genes = decode_str_array(self._h5["var"]["_index"][:])
        self._obs_cache: dict[str, ObsColumn] = {}
        self.matrix = (CsrSource if reader == "h5py" else AnndataSource)(path, layer)

    def has_obs(self, name: str) -> bool:
        return name in self._h5["obs"]

    def col(self, name: str) -> ObsColumn:
        if name not in self._obs_cache:
            if not self.has_obs(name):
                raise KeyError(
                    f"{self.name}: obs has no column {name!r}; "
                    f"available: {sorted(self._h5['obs'].keys())}"
                )
            self._obs_cache[name] = ObsColumn(self._h5["obs"][name])
        return self._obs_cache[name]

    def close(self) -> None:
        self.matrix.close()
        self._h5.close()


# --------------------------------------------------------------------------- #
# Selection plan
# --------------------------------------------------------------------------- #

class Block:
    """Cells of one group belonging to one output region."""

    __slots__ = ("dataset", "group", "condition", "sex", "rows", "barcodes")

    def __init__(self, dataset, group, condition, sex, rows, barcodes):
        self.dataset = dataset
        self.group = group
        self.condition = condition
        self.sex = sex
        self.rows = rows
        self.barcodes = barcodes

    def __len__(self) -> int:
        return len(self.rows)


def build_plan(datasets: dict[str, Dataset], specs, args):
    """region -> [Block, ...] plus a de-duplication report."""
    used: set[str] = set()
    plan: dict[str, list[Block]] = {}
    dropped: list[tuple[str, str, int]] = []

    for spec in specs:
        ds = datasets[spec["dataset"]]
        mask = (
            ds.col(args.status_column).match([spec["status"]])
            & ds.col(args.sex_column).match([spec["sex"]])
            & ds.col(args.region_column).match(spec["regions"])
        )
        pos = np.flatnonzero(mask)
        group = f"{spec['label']}_{spec['sex']}"
        if pos.size == 0:
            log(f"  {group}: no cells in {ds.name} for {list(spec['regions'])}")
            continue

        # drop cells already taken by an earlier group (shared controls)
        barcodes = ds.index[pos]
        keep = np.fromiter((b not in used for b in barcodes), dtype=bool, count=barcodes.size)
        n_dup = int((~keep).sum())
        if n_dup:
            dropped.append((group, ds.name, n_dup))
        pos, barcodes = pos[keep], barcodes[keep]
        used.update(barcodes.tolist())

        region_values = ds.col(args.region_column).take(pos)
        for region in spec["regions"]:
            sel = region_values == region
            if not sel.any():
                continue
            rows, bcs = pos[sel], barcodes[sel]          # ascending -> fast slicing
            if args.max_cells_per_group and rows.size > args.max_cells_per_group:
                rng = np.random.default_rng(args.seed)
                pick = np.sort(
                    rng.choice(rows.size, size=args.max_cells_per_group, replace=False)
                )
                rows, bcs = rows[pick], bcs[pick]
            plan.setdefault(region, []).append(
                Block(ds.name, group, spec["label"], spec["sex"], rows, bcs)
            )
    return plan, dropped


def sort_regions(regions) -> list[str]:
    order = {r: i for i, r in enumerate(REGION_ORDER)}
    return sorted(regions, key=lambda r: (order.get(r, len(order)), r))


def summarise_plan(plan, specs) -> pd.DataFrame:
    """region x group -> n_cells table."""
    groups = list(dict.fromkeys(f"{s['label']}_{s['sex']}" for s in specs))
    rows = []
    for region in sort_regions(plan.keys()):
        rec = {"region": region}
        total = 0
        for blk in plan[region]:
            rec[blk.group] = rec.get(blk.group, 0) + len(blk)
            total += len(blk)
        rec["total"] = total
        rows.append(rec)
    df = pd.DataFrame(rows)
    cols = ["region"] + [g for g in groups if g in df.columns] + ["total"]
    df = df[cols].fillna(0)
    for c in cols[1:]:
        df[c] = df[c].astype(int)
    return df


# --------------------------------------------------------------------------- #
# Loom writer
# --------------------------------------------------------------------------- #

def vstack_csr(mats) -> sp.csr_matrix:
    """Vertical concatenation of CSR blocks into one pre-allocated matrix.

    Avoids ``scipy.sparse.vstack`` (which goes through COO / bmat and needs
    several extra copies of a multi-hundred-million-nnz matrix).
    """
    n_rows = sum(m.shape[0] for m in mats)
    n_cols = mats[0].shape[1]
    nnz = sum(m.nnz for m in mats)
    big = max(nnz, n_rows, n_cols) >= np.iinfo(np.int32).max
    idx_dtype = np.int64 if big else np.int32

    data = np.empty(nnz, dtype=np.float32)
    indices = np.empty(nnz, dtype=idx_dtype)
    indptr = np.empty(n_rows + 1, dtype=idx_dtype)
    indptr[0] = 0

    row_off = nnz_off = 0
    for m in mats:
        if m.nnz:
            data[nnz_off:nnz_off + m.nnz] = m.data
            indices[nnz_off:nnz_off + m.nnz] = m.indices
        indptr[row_off + 1:row_off + 1 + m.shape[0]] = m.indptr[1:] + nnz_off
        row_off += m.shape[0]
        nnz_off += m.nnz

    return sp.csr_matrix((data, indices, indptr), shape=(n_rows, n_cols))


def count_detected_genes(matrix: sp.csr_matrix) -> np.ndarray:
    """Number of detected (>0) genes per cell, without a ``matrix > 0`` copy."""
    counts = matrix.getnnz(axis=1).astype(np.int32)
    if matrix.nnz and np.count_nonzero(matrix.data == 0):
        # rare: explicit zeros stored in the h5ad -> count them out properly
        rows = np.repeat(np.arange(matrix.shape[0], dtype=np.int64),
                         np.diff(matrix.indptr).astype(np.int64))
        counts = np.bincount(rows[matrix.data > 0], minlength=matrix.shape[0]).astype(np.int32)
    return counts


def _write_loom_attr(group, name, values, vlen_dtype) -> None:
    arr = np.asarray(values)
    if arr.dtype == object or arr.dtype.kind in ("U", "S"):
        group.create_dataset(name, data=np.array([str(x) for x in arr], dtype=object),
                             dtype=vlen_dtype)
    else:
        group.create_dataset(name, data=arr)


def create_loom_fast(path, matrix_csc, row_attrs, col_attrs, file_attrs=None,
                     chunk=1024, compression="gzip", compression_opts=2) -> None:
    """Write a Loom-3.0 file with h5py instead of ``loompy.create``.

    ``loompy.create`` stores a sparse matrix by densifying column windows into
    (64, 64) gzip chunks.  For a region-scale matrix (16k genes x 190k cells,
    ~450 M non-zeros) that is *millions* of tiny HDF5 chunk writes and takes
    hours.  The same data written with (1024, 1024) chunks takes a couple of
    minutes and is read back by ``loompy.connect`` / pySCENIC unchanged.
    """
    n_genes, n_cells = matrix_csc.shape
    vlen = h5py.special_dtype(vlen=str)

    with h5py.File(path, "w") as f:
        ds = f.create_dataset(
            "/matrix", shape=(n_genes, n_cells), dtype="float32",
            chunks=(min(chunk, n_genes), min(chunk, n_cells)),
            compression=compression, compression_opts=compression_opts, shuffle=True,
        )
        csc = matrix_csc.tocsc()
        window = max(1, min(4096, int(2e8 // max(1, n_genes))))
        for start in range(0, n_cells, window):
            stop = min(start + window, n_cells)
            ds[:, start:stop] = csc[:, start:stop].toarray()

        ra = f.create_group("/row_attrs")
        for name, values in row_attrs.items():
            _write_loom_attr(ra, name, values, vlen)
        ca = f.create_group("/col_attrs")
        for name, values in col_attrs.items():
            _write_loom_attr(ca, name, values, vlen)

        ga = f.create_group("/attrs")
        # loompy stores global attributes as scalars, so do the same
        ga["LOOM_SPEC_VERSION"] = "3.0.0"
        ga["CreationDate"] = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
        for key, value in (file_attrs or {}).items():
            ga[key] = str(value)


def write_region_loom(region, blocks, datasets, genes, out_path, args, build=True):
    """Concatenate the blocks of one region and write a loom (+ optional h5ad).

    With ``build=False`` nothing is written: the existing loom (if any) is only
    inspected so that ``--manifest-only`` can report / verify the outputs.
    """
    n_cells = sum(len(b) for b in blocks)
    groups = ";".join(sorted({b.group for b in blocks}))
    exists = os.path.exists(out_path)

    if not build or (exists and not args.overwrite):
        status, n_cells_file, n_genes_file = "missing", n_cells, len(genes)
        if exists:
            import loompy as lp

            with lp.connect(out_path, mode="r") as ds:
                n_genes_file, n_cells_file = int(ds.shape[0]), int(ds.shape[1])
            status = "ok" if n_cells_file == n_cells else "cell_count_mismatch"
        log(f"  {region}: {status:<20s} {out_path}")
        return dict(
            region=args.rename.get(region, region), region_source=region,
            file=out_path, status=status, n_cells=n_cells_file, n_genes=n_genes_file,
            size_gb=round(os.path.getsize(out_path) / 1e9, 3) if exists else 0.0,
            groups=groups,
        )

    log(f"  {region}: {n_cells:,} cells from {len(blocks)} block(s) -> {out_path}")

    matrices = []
    for blk in blocks:
        t0 = time.time()
        mat = datasets[blk.dataset].matrix.rows(blk.rows)
        log(f"    {blk.group:<8s} {len(blk):>8,} cells  {mat.nnz / 1e6:>8.1f} M nnz  "
            f"({time.time() - t0:.1f}s)")
        matrices.append(mat)

    t0 = time.time()
    matrix = matrices[0] if len(matrices) == 1 else vstack_csr(matrices)
    del matrices
    matrix = matrix.astype(np.float32, copy=False)
    log(f"    stacked: {matrix.shape[0]:,} x {matrix.shape[1]:,}  "
        f"{matrix.nnz / 1e6:.1f} M nnz ({time.time() - t0:.1f}s)")

    # ---- column attributes ------------------------------------------------
    t0 = time.time()
    cell_id = np.concatenate([b.barcodes for b in blocks]).astype(object)
    if len(np.unique(cell_id)) != len(cell_id):
        raise RuntimeError(f"{region}: duplicated barcodes in the loom")

    col_attrs = {
        "CellID": cell_id,
        "Group": np.concatenate([np.repeat(b.group, len(b)) for b in blocks]).astype(object),
        "Condition": np.concatenate([np.repeat(b.condition, len(b)) for b in blocks]).astype(object),
        "Sex": np.concatenate([np.repeat(b.sex, len(b)) for b in blocks]).astype(object),
        "Dataset": np.concatenate([np.repeat(b.dataset, len(b)) for b in blocks]).astype(object),
        "Region": np.repeat(args.rename.get(region, region), n_cells).astype(object),
        "RegionSource": np.repeat(region, n_cells).astype(object),
        "nGene": count_detected_genes(matrix),
        "nUMI": np.asarray(matrix.sum(axis=1)).ravel().astype(np.float32),
    }
    for col in args.celltype_columns:
        if not all(datasets[b.dataset].has_obs(col) for b in blocks):
            log(f"    [warn] obs column {col!r} missing in some input -- skipped")
            continue
        col_attrs[col] = np.concatenate(
            [datasets[b.dataset].col(col).take(b.rows) for b in blocks]
        ).astype(object)
    for attr, col in (("Sample", args.sample_column), ("Donor", args.donor_column)):
        if col and all(datasets[b.dataset].has_obs(col) for b in blocks):
            col_attrs[attr] = np.concatenate(
                [datasets[b.dataset].col(col).take(b.rows) for b in blocks]
            ).astype(object)
    log(f"    col_attrs built ({time.time() - t0:.1f}s)")

    # ---- write ------------------------------------------------------------
    tmp_path = out_path + ".tmp"
    if os.path.exists(tmp_path):
        os.remove(tmp_path)
    t0 = time.time()
    try:
        create_loom_fast(
            tmp_path,
            matrix.T,
            {"Gene": np.asarray(genes, dtype=object)},
            col_attrs,
            file_attrs={
                "note": "region-level subset loom (11subsetloom_meth.py)",
                "region": args.rename.get(region, region),
                "groups": ",".join(sorted({b.group for b in blocks})),
                "n_cells": str(n_cells),
                "matrix": args.layer or "X",
            },
        )
        os.replace(tmp_path, out_path)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

    size_gb = os.path.getsize(out_path) / 1e9
    log(f"    wrote {out_path}  ({size_gb:.2f} GB, {time.time() - t0:.1f}s)")

    # ---- optional h5ad for scanpy / R ------------------------------------
    if args.adata_output:
        import anndata as ad

        os.makedirs(args.adata_output, exist_ok=True)
        adata = ad.AnnData(
            X=matrix,
            obs=pd.DataFrame({k: np.asarray(v) for k, v in col_attrs.items()}, index=cell_id),
            var=pd.DataFrame(index=pd.Index(np.asarray(genes))),
        )
        adata_path = os.path.join(args.adata_output, f"{args.rename.get(region, region)}.h5ad")
        adata.write_h5ad(adata_path)
        log(f"    wrote {adata_path}")

    return dict(
        region=args.rename.get(region, region),
        region_source=region,
        file=out_path,
        status="written",
        n_cells=n_cells,
        n_genes=len(genes),
        size_gb=round(size_gb, 3),
        groups=";".join(sorted({b.group for b in blocks})),
    )


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="One loom per brain region for selected CUMS/CSDS/control groups.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input", action="append", metavar="NAME=H5AD", default=None,
                   help="dataset h5ad, repeatable (default: CUMS_4VN + CSDS_4VN)")
    p.add_argument("--output", default=DEFAULT_OUTPUT,
                   help=f"output directory for the region looms (default: {DEFAULT_OUTPUT})")
    p.add_argument("--adata-output", default=None,
                   help="optional directory for a per-region .h5ad copy")
    p.add_argument("--layer", default="counts",
                   help="h5ad layer to write ('counts' = raw UMIs, 'X' = normalised); "
                        "default: counts")
    p.add_argument("--config", default=None,
                   help="JSON file with optional keys inputs / groups / rename")
    p.add_argument("--groups", default=None,
                   help="comma separated subset of group keys, e.g. CUMS_M,CON_M")
    p.add_argument("--regions", default=None,
                   help="comma separated subset of output regions, e.g. AMY,STR")
    p.add_argument("--max_cells_per_group", type=int, default=0,
                   help="cap cells per group per region (0 = keep all)")
    p.add_argument("--seed", type=int, default=123, help="seed for the subsampling")
    p.add_argument("--status_column", default=STATUS_COL)
    p.add_argument("--sex_column", default=SEX_COL)
    p.add_argument("--region_column", default=REGION_COL)
    p.add_argument("--sample_column", default=SAMPLE_COL)
    p.add_argument("--donor_column", default=DONOR_COL)
    p.add_argument("--celltype_columns", default=",".join(CELLTYPE_COLS),
                   help="obs columns copied into col_attrs")
    p.add_argument("--rename", default=",".join(f"{k}={v}" for k, v in DEFAULT_RENAME.items()),
                   help="region rename for the output names, e.g. HPF=HIP")
    p.add_argument("--reader", choices=("h5py", "anndata"), default="h5py",
                   help="h5py = fast direct CSR reader (default); "
                        "anndata = backed reader for non-CSR matrices")
    p.add_argument("--overwrite", action="store_true", help="rebuild existing looms")
    p.add_argument("--dry-run", action="store_true", help="only print the cell-count table")
    p.add_argument("--manifest-only", action="store_true",
                   help="do not build anything: inspect the existing looms and rewrite "
                        "loom_manifest.csv / run_config.json")
    args = p.parse_args(argv)

    args.celltype_columns = [c.strip() for c in args.celltype_columns.split(",") if c.strip()]
    args.rename = dict(kv.split("=", 1) for kv in args.rename.split(",") if "=" in kv)

    if args.config:
        with open(args.config) as fh:
            cfg = json.load(fh)
        if "inputs" in cfg:
            args.input = [f"{k}={v}" for k, v in cfg["inputs"].items()]
        if "groups" in cfg:
            globals()["GROUP_SPECS"] = cfg["groups"]
        if "rename" in cfg:
            args.rename = cfg["rename"]

    return args


def resolve_inputs(spec_list):
    inputs = {}
    for spec in spec_list:
        if "=" in spec:
            name, path = spec.split("=", 1)
        else:
            name, path = os.path.splitext(os.path.basename(spec))[0], spec
        inputs[name] = path
    return inputs


def main(argv=None) -> int:
    args = parse_args(argv)

    # `kill -USR1 <pid>` dumps the Python stack -- handy when a big region seems
    # to stall, and harmless otherwise.
    try:
        import faulthandler
        import signal

        faulthandler.register(signal.SIGUSR1)
    except Exception:
        pass

    inputs = resolve_inputs(args.input) if args.input else dict(DEFAULT_INPUTS)

    specs = GROUP_SPECS
    if args.groups:
        wanted = {g.strip() for g in args.groups.split(",") if g.strip()}
        specs = [s for s in specs if f"{s['label']}_{s['sex']}" in wanted]
        missing = wanted - {f"{s['label']}_{s['sex']}" for s in specs}
        if missing:
            log(f"[warn] unknown group key(s): {sorted(missing)}")
    specs = [s for s in specs if s["dataset"] in inputs]
    if not specs:
        log("[error] no group left to build")
        return 2

    wanted_regions = None
    if args.regions:
        wanted_regions = {r.strip() for r in args.regions.split(",") if r.strip()}

    os.makedirs(args.output, exist_ok=True)
    group_keys = sorted({f"{s['label']}_{s['sex']}" for s in specs})

    log("=" * 78)
    log("region-level subset loom")
    log(f"  inputs : {', '.join(f'{k}={v}' for k, v in inputs.items())}")
    log(f"  output : {args.output}")
    log(f"  layer  : {args.layer}")
    log(f"  groups : {', '.join(group_keys)}")
    log(f"  cap    : {args.max_cells_per_group or 'no subsampling'} cells/group/region")
    log("=" * 78)

    datasets: dict[str, Dataset] = {}
    try:
        for name, path in inputs.items():
            if not any(s["dataset"] == name for s in specs):
                continue
            log(f"opening {name}: {path}")
            t0 = time.time()
            ds = Dataset(name, path, None if args.layer == "X" else args.layer, args.reader)
            datasets[name] = ds
            log(f"  {ds.index.size:,} cells x {ds.genes.size:,} genes "
                f"({time.time() - t0:.1f}s)")

        gene_sets = {tuple(ds.genes.tolist()) for ds in datasets.values()}
        if len(gene_sets) > 1:
            log("[warn] inputs have different var_names; the first file's genes are used")

        t0 = time.time()
        plan, dropped = build_plan(datasets, specs, args)
        for group, dsname, n in dropped:
            log(f"  {group}: dropped {n:,} already-used (control) cells from {dsname}")

        if not plan:
            log("[error] the selection is empty -- check the obs columns / values")
            return 2

        table = summarise_plan(plan, specs)
        if wanted_regions:
            table = table[table["region"].isin(wanted_regions)]
            plan = {r: b for r, b in plan.items() if args.rename.get(r, r) in wanted_regions}
            if not plan:
                log(f"[error] none of the requested regions {sorted(wanted_regions)} was found")
                return 2

        table = table.rename(columns={"region": "region(data)"})
        table.insert(1, "loom", [args.rename.get(r, r) for r in table["region(data)"]])
        log("")
        log("cells per region / group:")
        print(table.to_string(index=False), flush=True)
        log(f"total cells: {int(table['total'].sum()):,}  "
            f"(plan built in {time.time() - t0:.1f}s)")

        if args.dry_run:
            log("dry run -- nothing written")
            return 0

        genes = datasets[next(iter(datasets))].genes
        build = not args.manifest_only
        if not build:
            log("")
            log("manifest only -- existing looms are inspected, nothing is built")
        manifest = []
        for region in sort_regions(plan.keys()):
            out_path = os.path.join(args.output, f"{args.rename.get(region, region)}.loom")
            try:
                rec = write_region_loom(
                    region, plan[region], datasets, genes, out_path, args, build=build
                )
                if rec:
                    manifest.append(rec)
            except Exception as exc:                      # keep going for the other regions
                log(f"  [error] {region}: {type(exc).__name__}: {exc}")
                continue

        if manifest:
            mf = pd.DataFrame(manifest)[
                ["region", "region_source", "file", "status", "n_cells",
                 "n_genes", "size_gb", "groups"]
            ]
            mf_path = os.path.join(args.output, "loom_manifest.csv")
            mf.to_csv(mf_path, index=False)
            log("")
            log(mf.to_string(index=False))
            log(f"manifest: {mf_path}")

        with open(os.path.join(args.output, "run_config.json"), "w") as fh:
            json.dump(
                dict(
                    inputs=inputs, output=args.output, layer=args.layer,
                    groups=[{k: list(v) if isinstance(v, tuple) else v for k, v in s.items()}
                            for s in specs],
                    rename=args.rename, max_cells_per_group=args.max_cells_per_group,
                    seed=args.seed, celltype_columns=args.celltype_columns,
                ),
                fh, indent=2,
            )
        log("done")
        return 0
    finally:
        for ds in datasets.values():
            ds.close()


if __name__ == "__main__":
    sys.exit(main())
