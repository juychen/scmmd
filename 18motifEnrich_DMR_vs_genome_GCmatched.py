#!/usr/bin/env python
"""GC-matched genome-wide background: subsample genome windows so that their GC
distribution matches the DMR set bin-by-bin, then presence-based Fisher per motif.
Avoids the compositional (Simpson) artefacts of post-hoc GC adjustment.
"""
import numpy as np, pandas as pd
from scipy.stats import fisher_exact

W = '/tmp/enrich'
BASE = '/data2st2/junyi/output/atac1112/tobiasbam/BULK/BIND_DMR_enriched'
OUT = f'{BASE}/motif_enrichment_DMR_vs_genome_GCmatched.tsv'


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


dmr_len = read_fai(f'{W}/dmr8k.fa.fai')
gbg_len = read_fai(f'{W}/gbg.fa.fai')
dg = gc_of_fa(f'{W}/dmr8k.fa', dmr_len)
gg = gc_of_fa(f'{W}/gbg.fa', gbg_len)
bd = np.floor(dg * 100).astype(int)
bg = np.floor(gg * 100).astype(int)

rng = np.random.default_rng(20261112)
kd, kg = [], []
for b in sorted(set(bd) | set(bg)):
    D = list(dg.index[bd == b]); G = list(gg.index[bg == b])
    n = min(len(D), len(G))
    if n == 0:
        continue
    kd += list(rng.choice(D, n, replace=False))
    kg += list(rng.choice(G, n, replace=False))
print('matched windows: dmr', len(kd), 'genome', len(kg), flush=True)

dmr = load_hits(f'{W}/dmr8k_gc.tsv')
gbg = load_hits(f'{W}/gbg_gc.tsv')
kd = set(kd); kg = set(kg)
dmr = dmr[dmr.seq.isin(kd)]
gbg = gbg[gbg.seq.isin(kg)]
print('hits after matching:', len(dmr), len(gbg), flush=True)

dpres = dmr.drop_duplicates(['seq', 'motif']).groupby('motif')['seq'].agg(set)
gpres = gbg.drop_duplicates(['seq', 'motif']).groupby('motif')['seq'].agg(set)
n1, n2 = len(kd), len(kg)
rows = []
for m in sorted(set(dpres.index) | set(gpres.index)):
    a = len(dpres.get(m, set())); c = len(gpres.get(m, set()))
    orr, p = fisher_exact([[a, n1 - a], [c, n2 - c]], alternative='two-sided')
    rows.append((m, a, c, orr, p))
r = pd.DataFrame(rows, columns=['motif', 'dmr_pres', 'gbg_pres', 'or_fisher', 'p_fisher'])
r['pres_frac_dmr'] = r.dmr_pres / n1
r['pres_frac_gbg'] = r.gbg_pres / n2
r['fdr_fisher'] = bh(r.p_fisher)

bd_ = []
for reg in ['AMY', 'HIP', 'PFC']:
    d = pd.read_csv(f'{BASE}/{reg}_MC_footprints.bw/bindetect_results.txt', sep='\t')
    d = d[['name', 'MC_bound', 'total_tfbs', 'MC_MW_change', 'MC_MW_highlighted']].rename(columns={'name': 'motif'})
    d[f'{reg}_bf'] = d.MC_bound / d.total_tfbs.replace(0, np.nan)
    d = d.rename(columns={c: f'{reg}_{c}' for c in ['MC_bound', 'total_tfbs', 'MC_MW_change', 'MC_MW_highlighted']})
    bd_.append(d)
bd_ = bd_[0].merge(bd_[1], on='motif', how='outer').merge(bd_[2], on='motif', how='outer')
bd_['mean_bf'] = bd_[[f'{x}_bf' for x in ['AMY', 'HIP', 'PFC']]].mean(axis=1)
bd_['n_highlighted'] = bd_[[f'{x}_MC_MW_highlighted' for x in ['AMY', 'HIP', 'PFC']]].sum(axis=1)

tf = pd.read_csv(f'{BASE}/motif_tf_annotation.tsv', sep='\t')
res = r.merge(bd_, on='motif', how='left').merge(tf, on='motif', how='left')
res = res.sort_values('p_fisher')
res.to_csv(OUT, sep='\t', index=False)
print('wrote', OUT, res.shape, flush=True)

for a in [0.05, 0.01]:
    print(f'FDR<{a} OR>1:', ((res.fdr_fisher < a) & (res.or_fisher > 1)).sum(),
          ' OR<1:', ((res.fdr_fisher < a) & (res.or_fisher < 1)).sum(), flush=True)

cols = ['motif', 'tf', 'pres_frac_dmr', 'pres_frac_gbg', 'or_fisher', 'fdr_fisher',
        'AMY_bf', 'HIP_bf', 'PFC_bf', 'mean_bf', 'n_highlighted']
print(res[(res.fdr_fisher < 0.05) & (res.or_fisher > 1)].nlargest(20, 'or_fisher')[cols].to_string(index=False))
