import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

B = '/data2st2/junyi/output/atac1112/tobiasbam/BULK/BIND_DMR_enriched/'
m = pd.read_csv(B + 'motif_enrichment_DMR_vs_genome_GCmatched.tsv', sep='\t')
f = pd.read_csv(B + 'motif_enrichment_DMR_flank_genome.tsv', sep='\t')
t = pd.read_csv(B + 'tf_enrichment_summary.tsv', sep='\t')
j = m.merge(f[['motif', 'or_mh_flank', 'fdr_mh_flank']], on='motif')

en_g = (j.fdr_fisher < 0.05) & (j.or_fisher > 1)
en_f = (j.fdr_mh_flank < 0.05) & (j.or_mh_flank > 1)
both = en_g & en_f
gonly = en_g & ~en_f
fonly = en_f & ~en_g

fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.8))
ax = axes[0]
x = np.log2(j.or_mh_flank.clip(.2, 8)); y = np.log2(j.or_fisher.clip(.2, 8))
ax.scatter(x[~en_g & ~en_f], y[~en_g & ~en_f], s=9, c='#cfcfcf', alpha=.45, lw=0, label='not enriched')
ax.scatter(x[fonly], y[fonly], s=16, c='#3b7dd8', alpha=.8, lw=0, label='flank only (n=%d)' % fonly.sum())
ax.scatter(x[gonly], y[gonly], s=16, c='#e8a33d', alpha=.8, lw=0, label='genome only (n=%d)' % gonly.sum())
ax.scatter(x[both], y[both], s=16, c='#d62728', alpha=.65, lw=0, label='both (n=%d)' % both.sum())
lim = [-2.4, 2.2]
ax.plot(lim, lim, ls='--', c='k', lw=.7)
ax.axhline(0, c='k', lw=.4); ax.axvline(0, c='k', lw=.4)
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel('log2 OR  (DMR vs local flanks)')
ax.set_ylabel('log2 OR  (DMR vs GC-matched genome)')
ax.set_title('Motif enrichment: two backgrounds')
labs = {
    'taipale_cyt_meth_GMEB1_NKACGTNNNNACGTMN_FL_meth_repr': ('Gmeb1', (6, 6)),
    'transfac_public_M00039': ('Creb1', (-6, 11)),
    'jaspar_MA0604.1': ('Atf1', (7, -4)),
    'homer_NANTGCSGCA_Zfp57': ('Zfp57', (6, -3)),
    'metacluster_67.2': ('E2f2/3 (mc67.2)', (-125, 34)),
    'metacluster_55.3': ('E2f6/Nrf1 (mc55.3)', (-125, 10)),
    'metacluster_161.40': ('Msx2', (4, 6)),
    'swissregulon_hs_PAX2': ('Pax2', (7, 2)),
    'metacluster_71.5': ('Hif1a', (-30, 4)),
    'cisbp_M07803': ('Myc', (-18, -12)),
    'hocomoco_ZBT14_HUMAN.H11MO.0.C': ('Zbtb14', (10, -22)),
    'metacluster_131.6': ('Egr2/3/4', (-50, -8)),
    'c2h2_zfs_M0393': ('Klf9', (-26, 6)),
    'taipale_tf_pairs_E2F3_HES7_NNNNGCGCSNNNNNCACGTGNN_CAP': ('E2f3/Hes7', (8, 4)),
}
for mo, (lab, off) in labs.items():
    r = j[j.motif == mo]
    if len(r) == 0:
        continue
    r = r.iloc[0]
    kw = dict(fontsize=8, xytext=off, textcoords='offset points', color='#111')
    if lab.startswith('E2f'):
        kw['arrowprops'] = dict(arrowstyle='-', lw=.6, color='#666')
    ax.annotate(lab, (np.log2(max(r.or_mh_flank, .2)), np.log2(max(r.or_fisher, .2))), **kw)
ax.legend(fontsize=8, frameon=False, loc='lower right')

ax = axes[1]
top = t.nlargest(15, 'n_enriched').iloc[::-1]
ax.barh(range(len(top)), top.n_enriched, color='#4c78a8', height=.72)
ax.set_yticks(range(len(top)))
ax.set_yticklabels(['%s (%d/%d)' % (r.tf, r.n_enriched, r.n_motifs) for r in top.itertuples()], fontsize=8.5)
ax.set_xlabel('# enriched motifs (FDR<0.05, OR>1)')
ax.set_title('TFs with most enriched motifs (genome background)')
for i, r in enumerate(top.itertuples()):
    ax.text(r.n_enriched + .3, i, 'max OR %.2f' % r.best_or, va='center', fontsize=7.2, color='#444')
ax.set_xlim(0, top.n_enriched.max() * 1.22)

plt.tight_layout()
plt.savefig('/tmp/enrich/enrich_genome_TF.png', dpi=200)
print('saved figure')

strong = j[both & (j.mean_bf > 0.10)].copy().sort_values('mean_bf', ascending=False)
out = strong[['motif', 'tf', 'pres_frac_dmr', 'pres_frac_gbg', 'or_fisher', 'fdr_fisher',
              'or_mh_flank', 'fdr_mh_flank', 'AMY_bf', 'HIP_bf', 'PFC_bf', 'mean_bf',
              'n_highlighted', 'AMY_MC_MW_change', 'HIP_MC_MW_change', 'PFC_MC_MW_change']]
out.to_csv(B + 'top_motif_evidence.tsv', sep='\t', index=False)
print('strong table', out.shape)
