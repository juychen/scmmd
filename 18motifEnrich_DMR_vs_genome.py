#!/usr/bin/env python
"""Motif enrichment in DMRs against two backgrounds:
  (1) local 500bp flanks (previous run), (2) genome-wide random 500bp windows.
Presence-based Fisher + GC-stratified Cochran-Mantel-Haenszel, BH FDR; merged with
BINDetect binding metrics and the cisTarget motif->TF annotation.
"""
import numpy as np, pandas as pd
from scipy.stats import fisher_exact, chi2

W = '/tmp/enrich'
BASE = '/data2st2/junyi/output/atac1112/tobiasbam/BULK/BIND_DMR_enriched'
OUT = f'{BASE}/motif_enrichment_DMR_flank_genome.tsv'
TF = f'{BASE}/motif_tf_annotation.tsv'


def load_hits(path):
    return pd.read_csv(path, sep='\t', header=None, usecols=[0, 3],
                       names=['seq', 'motif'], dtype=str, engine='c')


def read_fai(p):
    return {ln.split('\t')[0]: int(ln.split('\t')[1]) for ln in open(p)}


def gc_of_fa(path, keep):
    gc = {}; name = None; s = 0; n = 0
    for line in open(path):
        if line[0] == '>':
            if name is not None:
                gc[name] = s / n if n else np.nan
            name = line[1:].strip().split()[0]; s = 0; n = 0
        else:
            up = line.strip().upper(); s += up.count('G') + up.count('C'); n += len(up)
    if name is not None:
        gc[name] = s / n if n else np.nan
    return pd.Series(gc).reindex(list(keep))


def bh(p):
    p = np.asarray(p, float); n = len(p); q = np.full(n, np.nan)
    ok = ~np.isnan(p); pp = p[ok]; o = np.argsort(pp)
    r = pp[o] * n / (np.arange(len(pp)) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    qq = np.empty(len(pp)); qq[o] = np.minimum(r, 1); q[ok] = qq
    return q


def enrich(dmr_hits, dmr_pres, dmr_len, bg_hits, bg_pres, bg_len, bg_name):
    """Return per-motif stats vs one background group."""
    n1, n2 = len(dmr_len), len(bg_len)
    dmr_gc = gc_of_fa(f'{W}/dmr8k.fa', dmr_len)
    bg_gc = gc_of_fa(f'{W}/{bg_name}.fa', bg_len)
    pool = pd.concat([dmr_gc, bg_gc])
    qs = np.quantile(pool, np.linspace(0, 1, 11)); qs[0] = -1; qs[-1] = 2
    ds = pd.cut(dmr_gc, qs, labels=False, include_lowest=True).astype(int)
    bs = pd.cut(bg_gc, qs, labels=False, include_lowest=True).astype(int)
    D = {s: set(ds[ds == s].index) for s in range(10) if (ds == s).any()}
    B = {s: set(bs[bs == s].index) for s in range(10) if (bs == s).any()}
    dset = dmr_pres.groupby('motif')['seq'].agg(set) if 'motif' in dmr_pres else dmr_pres
    bset = bg_pres.groupby('motif')['seq'].agg(set) if 'motif' in bg_pres else bg_pres
    motifs = sorted(set(dmr_hits.index) | set(bg_hits.index))
    rows = []
    for m in motifs:
        dh = dset.get(m, set()); bh_ = bset.get(m, set())
        a, c = len(dh), len(bh_)
        or_raw, p_raw = fisher_exact([[a, n1 - a], [c, n2 - c]], alternative='two-sided')
        num = den = oe = var = 0.0
        for s in range(10):
            if s not in D or s not in B:
                continue
            N1, N2 = len(D[s]), len(B[s])
            if N1 == 0 or N2 == 0:
                continue
            as_ = len(dh & D[s]); cs_ = len(bh_ & B[s]); bs_ = N1 - as_; ds_ = N2 - cs_
            N = N1 + N2
            num += as_ * ds_ / N; den += bs_ * cs_ / N
            E = N1 * (as_ + cs_) / N
            V = N1 * N2 * (as_ + cs_) * (bs_ + ds_) / (N * N * (N - 1)) if N > 1 else 0
            oe += as_ - E; var += V
        or_mh = num / den if den > 0 else np.nan
        p_mh = chi2.sf(oe ** 2 / var, 1) if var > 0 else np.nan
        rows.append((m, a, c, dmr_hits.get(m, 0), bg_hits.get(m, 0), or_raw, p_raw, or_mh, p_mh))
    r = pd.DataFrame(rows, columns=['motif', f'dmr_pres_{bg_name}', f'{bg_name}_pres',
                                    f'dmr_hits_{bg_name}', f'{bg_name}_hits',
                                    f'or_raw_{bg_name}', f'p_raw_{bg_name}',
                                    f'or_mh_{bg_name}', f'p_mh_{bg_name}'])
    r[f'pres_frac_dmr_{bg_name}'] = r[f'dmr_pres_{bg_name}'] / n1
    r[f'pres_frac_{bg_name}'] = r[f'{bg_name}_pres'] / n2
    r[f'fdr_raw_{bg_name}'] = bh(r[f'p_raw_{bg_name}'])
    r[f'fdr_mh_{bg_name}'] = bh(r[f'p_mh_{bg_name}'])
    return r


print('loading hits ...', flush=True)
dmr = load_hits(f'{W}/dmr8k_gc.tsv')
flk = load_hits(f'{W}/flank_gc.tsv')
gbg = load_hits(f'{W}/gbg_gc.tsv')
print(len(dmr), len(flk), len(gbg), flush=True)

dmr_len = read_fai(f'{W}/dmr8k.fa.fai')
flk_len = read_fai(f'{W}/flank.fa.fai')
gbg_len = read_fai(f'{W}/gbg.fa.fai')
print('bp', sum(dmr_len.values()), sum(flk_len.values()), sum(gbg_len.values()), flush=True)

dmr_hits = dmr['motif'].value_counts()
dmr_pres = dmr.drop_duplicates(['seq', 'motif'])
flk_hits = flk['motif'].value_counts()
flk_pres = flk.drop_duplicates(['seq', 'motif'])
gbg_hits = gbg['motif'].value_counts()
gbg_pres = gbg.drop_duplicates(['seq', 'motif'])

res = enrich(dmr_hits, dmr_pres, dmr_len, flk_hits, flk_pres, flk_len, 'flank')
res2 = enrich(dmr_hits, dmr_pres, dmr_len, gbg_hits, gbg_pres, gbg_len, 'gbg')
res = res.drop(columns=[c for c in res.columns if c.endswith('_flank') and c.startswith('dmr_')]).merge(
    res2.drop(columns=[c for c in res2.columns if c.endswith('_gbg') and c.startswith('dmr_')]), on='motif')

# BINDetect binding metrics
bd = []
for reg in ['AMY', 'HIP', 'PFC']:
    d = pd.read_csv(f'{BASE}/{reg}_MC_footprints.bw/bindetect_results.txt', sep='\t')
    d = d[['name', 'total_tfbs', 'MC_mean_score', 'MC_bound', 'MW_bound', 'MC_MW_change',
           'MC_MW_pvalue', 'MC_MW_highlighted']].rename(columns={'name': 'motif'})
    d[f'{reg}_bf'] = d.MC_bound / d.total_tfbs.replace(0, np.nan)
    d = d.rename(columns={c: f'{reg}_{c}' for c in ['total_tfbs', 'MC_mean_score', 'MC_bound',
                                                    'MW_bound', 'MC_MW_change', 'MC_MW_pvalue',
                                                    'MC_MW_highlighted']})
    bd.append(d)
bd = bd[0].merge(bd[1], on='motif', how='outer').merge(bd[2], on='motif', how='outer')
bd['max_abs_change'] = bd[[f'{r}_MC_MW_change' for r in ['AMY', 'HIP', 'PFC']]].abs().max(axis=1)
bd['n_highlighted'] = bd[[f'{r}_MC_MW_highlighted' for r in ['AMY', 'HIP', 'PFC']]].sum(axis=1)

tf = pd.read_csv(TF, sep='\t') if __import__('os').path.exists(TF) else None
out = res.merge(bd, on='motif', how='left')
if tf is not None:
    out = out.merge(tf, on='motif', how='left')
out = out.sort_values('p_mh_gbg', na_position='last')
out.to_csv(OUT, sep='\t', index=False)
print('wrote', OUT, out.shape, flush=True)

for name in ['flank', 'gbg']:
    sig = out[(out[f'fdr_mh_{name}'] < 0.05) & (out[f'or_mh_{name}'] > 1)]
    dep = out[(out[f'fdr_mh_{name}'] < 0.05) & (out[f'or_mh_{name}'] < 1)]
    print(f'[{name}] enriched {len(sig)}, depleted {len(dep)}', flush=True)

cols = ['motif', 'tf', 'pres_frac_dmr_gbg', 'pres_frac_gbg', 'or_mh_gbg', 'fdr_mh_gbg',
        'or_mh_flank', 'fdr_mh_flank', 'AMY_bf', 'HIP_bf', 'PFC_bf', 'max_abs_change']
print(out.nlargest(15, 'or_mh_gbg')[cols].to_string(index=False), flush=True)
