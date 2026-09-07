#!/usr/bin/env python3
"""Compare SCENIC ctx regulons between L1 (coarse cell types) and L2 (fine subtypes).

For every L2 subtype we map it back to its parent L1 category, take the
significant regulons (max NES over motifs > NES_THRESHOLD), and measure how
consistent the active TF sets and NES profiles are between the parent L1 run
and the L2 subtypes run.
"""
import glob
import os
import re

import numpy as np
import pandas as pd
from scipy import stats

L1_DIR = "/data1st1/junyi/output/sn0827L1/scenic_ctx"
L2_DIR = "/data1st1/junyi/output/sn0827/scenic_ctx"
OUT_DIR = "/home/junyichen/code/scmmd/output/l1_l2_ctx_consistency"
NES_THRESHOLD = 3.0  # pySCENIC convention for significant regulons
TOP_K = 20

# subtype name -> L1 category (region prefix is stripped before matching here)
SUBTYPE2CATEGORY = {
    "Astrocyte": "Astrocyte",
    "Choroid_plexus_cell": "Epen",
    "Ependymal_cell": "Epen",
    "Hypendymal_cell": "Epen",
    "Tanycyte": "Epen",
    "Microglia": "Microglia",
    "MOL": "Oligo",
    "NFOL": "Oligo",
    "MFOL": "Oligo",
    "COP": "Oligo",
    "Immature_cell": "Oligo",
    "OPC": "OPC",
    "Endothelial_cell": "Vascular",
    "Pericyte": "Vascular",
    "VSMC": "Vascular",
    "VLMC": "Vascular",
    "Arachnoid_Barrier_cell": "Vascular",
    "Lymphocyte": "Immune",
    "Perivascular_Macrophage": "Immune",
    "Chol": "Chol",
    "Hist": "Hist",
    "Dopa": "Dopa",
    "Sero": "Sero",
}


def load_ctx(path):
    try:
        df = pd.read_csv(path, skiprows=3, header=None,
                         names=["TF", "MotifID", "AUC", "NES", "MotifSimilarityQvalue",
                                "OrthologousIdentity", "Annotation", "Context",
                                "TargetGenes", "RankAtMax"])
    except pd.errors.EmptyDataError:
        return pd.Series(dtype=float), set()
    # one TF can have several motifs; keep its best NES
    tf_nes = df.groupby("TF")["NES"].max()
    active = set(tf_nes[tf_nes > NES_THRESHOLD].index)
    return tf_nes, active


def parse_l1_name(fname):
    m = re.match(r"ctx_([A-Z]+)_(.+)\.csv$", fname)
    return (m.group(1), m.group(2)) if m else (None, None)


def subtype_to_category(subtype, region):
    # some subtype names carry the region prefix again, e.g. "HY_Hist"
    if subtype.startswith(region + "_"):
        subtype = subtype[len(region) + 1:]
    if subtype.endswith("_Glut"):
        return "Glut"
    if subtype.endswith("_GABA"):
        return "GABA"
    base = re.sub(r"-\d+$", "", subtype)
    if base in SUBTYPE2CATEGORY:
        return SUBTYPE2CATEGORY[base]
    if base in SUBTYPE2CATEGORY.values():  # already a category name
        return base
    return None


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    l1_files = sorted(glob.glob(os.path.join(L1_DIR, "ctx_*.csv")))
    l2_files = sorted(glob.glob(os.path.join(L2_DIR, "ctx_*.csv")))

    l1_data, l2_map = {}, {}
    for f in l2_files:
        fname = os.path.basename(f)
        region, subtype = parse_l1_name(fname)
        cat = subtype_to_category(subtype, region)
        if cat is None:
            print(f"[warn] unmapped L2 subtype: {fname}")
            continue
        l2_map.setdefault((region, cat), []).append((subtype, f))

    rows = []
    for f in l1_files:
        region, cat = parse_l1_name(os.path.basename(f))
        key = (region, cat)
        if key not in l2_map:
            print(f"[warn] no L2 subtypes matched for {key}")
            continue
        l1_nes, l1_active = load_ctx(f)

        subtype_active, subtype_nes = {}, {}
        for subtype, sf in l2_map[key]:
            nes, active = load_ctx(sf)
            subtype_active[subtype] = active
            subtype_nes[subtype] = nes

        union = set().union(*subtype_active.values())
        inter = set.intersection(*subtype_active.values()) if subtype_active else set()

        # mean NES profile over subtypes (union of TFs)
        l2_mean = pd.DataFrame({s: n for s, n in subtype_nes.items()}).mean(axis=1)

        jaccard_union = len(l1_active & union) / len(l1_active | union) if l1_active | union else np.nan
        jaccard_inter = len(l1_active & inter) / len(l1_active | inter) if l1_active | inter else np.nan

        top_l1 = l1_nes.sort_values(ascending=False).head(TOP_K).index

        def top_tfs(s):
            nes_s = subtype_nes[s]
            act = subtype_active[s]
            return set(sorted(act, key=lambda t: -nes_s[t])[:TOP_K])

        subtype_top = {s: top_tfs(s) for s in subtype_active}
        top_all = set().union(*subtype_top.values()) if subtype_top else set()
        recall_any = np.mean([tf in union for tf in top_l1])
        top_rank_recall = np.mean([tf in top_all for tf in top_l1])
        top_in_any = np.mean([any(tf in subtype_top[s] for s in subtype_top) for tf in top_l1])

        common = l2_mean.index.intersection(l1_nes.index)
        rho, pval = stats.spearmanr(l1_nes[common], l2_mean[common])

        rows.append({
            "region": region, "category": cat,
            "n_subtypes": len(subtype_active),
            "n_active_L1": len(l1_active),
            "n_active_L2_union": len(union),
            "n_active_L2_intersection": len(inter),
            "jaccard_L1_vs_L2union": jaccard_union,
            "jaccard_L1_vs_L2intersection": jaccard_inter,
            f"recall_L1_top{TOP_K}_in_L2union": recall_any,
            f"recall_L1_top{TOP_K}_in_L2top{TOP_K}": top_rank_recall,
            "frac_L1_top_topK_in_some_subtype": top_in_any,
            "spearman_nes_L1_vs_L2mean": rho,
            "spearman_p": pval,
        })

    res = pd.DataFrame(rows).sort_values(["region", "category"])
    res.to_csv(os.path.join(OUT_DIR, "consistency_per_category.csv"), index=False)

    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 20)
    print(res.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    print("\n=== overall ===")
    for col in ["jaccard_L1_vs_L2union", f"recall_L1_top{TOP_K}_in_L2union",
                f"recall_L1_top{TOP_K}_in_L2top{TOP_K}",
                "spearman_nes_L1_vs_L2mean"]:
        s = res[col].astype(float)
        print(f"{col}: median={s.median():.3f}  mean={s.mean():.3f}  min={s.min():.3f}  max={s.max():.3f}")


if __name__ == "__main__":
    main()
