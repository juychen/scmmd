#!/usr/bin/env python
"""Motif enrichment: DMR windows vs 500bp flanks (presence-based, GC-stratified).
Inputs (prepared in /tmp/enrich): dmr8k.tsv, flank.tsv (TOBIAS TFBScan output),
dmr8k.fa(.fai), flank.fa(.fai), plus BINDetect bindetect_results.txt of 3 regions.
"""
import numpy as np, pandas as pd
from scipy.stats import fisher_exact, chi2

W = '/tmp/enrich'
OUT = '/data2st2/junyi/output/atac1112/tobiasbam/BULK/BIND_DMR_enriched/motif_enrichment_DMR_vs_flank.tsv'


def load_hits(path):
    return pd.read_csv(path, sep='\t', header=None, usecols=[0, 3],
                       names=['seq', 'motif'], dtype=str, engine='c')


print('loading hits ...', flush=True)
dmr = load_hits(f'{W}/dmr8k.tsv')
flk = load_hits(f'{W}/flank.tsv')
print(len(dmr), len(flk), flush=True)

dmr_hits = dmr['motif'].value_counts()
flk_hits = flk['motif'].value_counts()
dmr_pairs = dmr.drop_duplicates(['seq', 'motif'])
flk_pairs = flk.drop_duplicates(['seq', 'motif'])
dmr_pres = dmr_pairs['motif'].value_counts()
flk_pres = flk_pairs['motif'].value_counts()


def read_fai(p):
    d = {}
    for line in open(p):
        f = line.rstrip('\n').split('\t')
        d[f[0]] = int(f[1])
    return d


dmr_len = read_fai(f'{W}/dmr8k.fa.fai')
flk_len = read_fai(f'{W}/flank.fa.fai')
dmr_bp, flk_bp = sum(dmr_len.values()), sum(flk_len.values())
print('seqs/bp', len(dmr_len), dmr_bp, len(flk_len), flk_bp, flush=True)


def gc_of_fa(path, keep):
    gc = {}; name = None; s = 0; n = 0
    with open(path) as fh:
        for line in fh:
            if line[0] == '>':
                if name is not None:
                    gc[name] = s / n if n else np.nan
                name = line[1:].strip().split()[0]; s = 0; n = 0
            else:
                up = line.strip().upper()
                s += up.count('G') + up.count('C'); n += len(up)
        if name is not None:
            gc[name] = s / n if n else np.nan
    return pd.Series(gc).reindex(list(keep))


dmr_gc = gc_of_fa(f'{W}/dmr8k.fa', dmr_len)
flk_gc = gc_of_fa(f'{W}/flank.fa', flk_len)

pool = pd.concat([dmr_gc, flk_gc])
qs = np.quantile(pool, np.linspace(0, 1, 11))
qs[0] = -1; qs[-1] = 2
dmr_s = pd.cut(dmr_gc, qs, labels=False, include_lowest=True).astype(int)
flk_s = pd.cut(flk_gc, qs, labels=False, include_lowest=True).astype(int)
D = {s: set(dmr_s[dmr_s == s].index) for s in range(10) if (dmr_s == s).any()}
F = {s: set(flk_s[flk_s == s].index) for s in range(10) if (flk_s == s).any()}

dmr_sets = dmr_pairs.groupby('motif')['seq'].agg(set)
flk_sets = flk_pairs.groupby('motif')['seq'].agg(set)

motifs = sorted(set(dmr_hits.index) | set(flk_hits.index))
print('motifs', len(motifs), flush=True)

n1, n2 = len(dmr_len), len(flk_len)
rows = []
for m in motifs:
    dh = dmr_sets.get(m, set()); fh = flk_sets.get(m, set())
    a, c = len(dh), len(fh)
    or_raw, p_raw = fisher_exact([[a, n1 - a], [c, n2 - c]], alternative='two-sided')
    num = den = 0.0; oe = 0.0; var = 0.0
    for s in range(10):
        if s not in D or s not in F:
            continue
        N1, N2 = len(D[s]), len(F[s])
        if N1 == 0 or N2 == 0:
            continue
        as_ = len(dh & D[s]); cs_ = len(fh & F[s])
        bs = N1 - as_; ds = N2 - cs_; N = N1 + N2
        num += as_ * ds / N; den += bs * cs_ / N
        E = N1 * (as_ + cs_) / N
        V = N1 * N2 * (as_ + cs_) * (bs + ds) / (N * N * (N - 1)) if N > 1 else 0
        oe += as_ - E; var += V
    or_mh = num / den if den > 0 else np.nan
    p_mh = chi2.sf(oe ** 2 / var, 1) if var > 0 else np.nan
    rows.append((m, a, c, dmr_hits.get(m, 0), flk_hits.get(m, 0), or_raw, p_raw, or_mh, p_mh))

res = pd.DataFrame(rows, columns=['motif', 'dmr_pres', 'flank_pres', 'dmr_hits', 'flank_hits',
                                  'or_raw', 'p_raw', 'or_mh', 'p_mh'])
res['pres_frac_dmr'] = res.dmr_pres / n1
res['pres_frac_flank'] = res.flank_pres / n2
res['hits_pkb_dmr'] = res.dmr_hits / dmr_bp * 1000
res['hits_pkb_flank'] = res.flank_hits / flk_bp * 1000
res['hit_rate_ratio'] = (res.dmr_hits / dmr_bp) / (res.flank_hits / flk_bp).replace(0, np.nan)


def bh(p):
    p = np.asarray(p, float); n = len(p); q = np.full(n, np.nan)
    ok = ~np.isnan(p); pp = p[ok]; o = np.argsort(pp)
    ranked = pp[o] * n / (np.arange(len(pp)) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    qq = np.empty(len(pp)); qq[o] = np.minimum(ranked, 1)
    q[ok] = qq
    return q


res['fdr_raw'] = bh(res.p_raw)
res['fdr_mh'] = bh(res.p_mh)

base = '/data2st2/junyi/output/atac1112/tobiasbam/BULK/BIND_DMR_enriched'
bd = []
for reg in ['AMY', 'HIP', 'PFC']:
    d = pd.read_csv(f'{base}/{reg}_MC_footprints.bw/bindetect_results.txt', sep='\t')
    d = d[['name', 'total_tfbs', 'MC_mean_score', 'MC_bound', 'MW_mean_score', 'MW_bound',
           'MC_MW_change', 'MC_MW_pvalue', 'MC_MW_highlighted']].rename(columns={'name': 'motif'})
    for c in ['total_tfbs', 'MC_mean_score', 'MC_bound', 'MW_mean_score', 'MW_bound',
              'MC_MW_change', 'MC_MW_pvalue', 'MC_MW_highlighted']:
        d[f'{reg}_{c}'] = d.pop(c)
    bd.append(d)
bd = bd[0].merge(bd[1], on='motif', how='outer').merge(bd[2], on='motif', how='outer')
for reg in ['AMY', 'HIP', 'PFC']:
    bd[f'{reg}_bf'] = bd[f'{reg}_MC_bound'] / bd[f'{reg}_total_tfbs'].replace(0, np.nan)
bd['max_abs_change'] = bd[[f'{r}_MC_MW_change' for r in ['AMY', 'HIP', 'PFC']]].abs().max(axis=1)
bd['n_highlighted'] = bd[[f'{r}_MC_MW_highlighted' for r in ['AMY', 'HIP', 'PFC']]].sum(axis=1)

out = res.merge(bd, on='motif', how='left').sort_values(['p_mh', 'p_raw'], na_position='last')
out.to_csv(OUT, sep='\t', index=False)
print('wrote', OUT, out.shape, flush=True)

print('enriched FDR<0.05:', int(((out.fdr_mh < 0.05) & (out.or_mh > 1)).sum()),
      ' depleted:', int(((out.fdr_mh < 0.05) & (out.or_mh < 1)).sum()), flush=True)
print('\ntop by MH OR (FDR<0.05):')
cols = ['motif', 'pres_frac_dmr', 'pres_frac_flank', 'or_mh', 'fdr_mh', 'hit_rate_ratio',
        'AMY_bf', 'HIP_bf', 'PFC_bf', 'max_abs_change']
print(out[out.fdr_mh < 0.05].nlargest(15, 'or_mh')[cols].to_string(index=False))
print('\ntop by raw OR:')
print(out.nlargest(10, 'or_raw')[['motif', 'pres_frac_dmr', 'pres_frac_flank', 'or_raw', 'fdr_raw']].to_string(index=False))
print('\nspearman or_mh vs bf:',
      out[['or_mh', 'AMY_bf', 'HIP_bf', 'PFC_bf']].corr(method='spearman').loc['or_mh'].to_dict(), flush=True)
print('spearman or_mh vs mean score:',
      out[['or_mh', 'AMY_MC_mean_score', 'HIP_MC_mean_score', 'PFC_MC_mean_score']].corr(method='spearman').loc['or_mh'].to_dict(), flush=True)
