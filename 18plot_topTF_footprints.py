import os, shutil, numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE = '/data2st2/junyi/output/atac1112/tobiasbam/BULK/BIND_DMR_enriched'
FIG = os.path.join(BASE, 'figures')
FP = os.path.join(FIG, 'footprints')
os.makedirs(FP, exist_ok=True)

# ---- 1. copy summary figures -------------------------------------------------
for src, dst in [
    ('/tmp/enrich/enrich_overview.png', '01_enrich_overview_flank.png'),
    ('/tmp/enrich/enrich_genome_TF.png', '02_enrich_two_backgrounds_TF.png'),
    ('/tmp/enrich/footprint_support.png', '03_footprint_support.png'),
]:
    if os.path.exists(src):
        shutil.copy2(src, os.path.join(FIG, dst))
        print('copied', dst)

# ---- 2. tables ---------------------------------------------------------------
m = pd.read_csv(os.path.join(BASE, 'motif_enrichment_DMR_vs_genome_GCmatched.tsv'), sep='\t')
ann = pd.read_csv(os.path.join(BASE, 'motif_tf_annotation.tsv'), sep='\t')
tfsum = pd.read_csv(os.path.join(BASE, 'tf_enrichment_summary.tsv'), sep='\t')

ex = ann[['motif']].assign(tf=ann.tf.str.split(';')).explode('tf')
j = ex.merge(m.drop(columns=['tf']), on='motif', how='left')
en = j[(j.fdr_fisher < 0.05) & (j.or_fisher > 1)]
ntf_map = dict(zip(ann.motif, ann.n_tf))
regions = ['AMY', 'HIP', 'PFC']

top = tfsum.nlargest(15, 'n_enriched').reset_index(drop=True)

def read_signal(region, motif):
    p = os.path.join(BASE, f'{region}_MC_footprints.bw', motif, 'beds', f'{motif}_all.bed_signal.txt')
    mc = mw = None
    for line in open(p):
        if line.startswith('#'):
            continue
        parts = line.rstrip('\n').split('\t')
        vals = np.array(parts[2].split(','), dtype=float)
        if '_MC_' in parts[0]:
            mc = vals
        elif '_MW_' in parts[0]:
            mw = vals
    return mc, mw

used, rows = set(), []
for _, r in top.iterrows():
    name = r['tf']
    cand = en[en.tf == name].sort_values(['or_fisher', 'mean_bf'], ascending=False)
    cand = cand[~cand.motif.isin(used)]
    if len(cand) == 0:
        cand = m[(m.motif.isin(ann[ann.tf.str.split(';').apply(lambda x: name in x)].motif))].sort_values('mean_bf', ascending=False)
    if len(cand) == 0:
        print('skip', name); continue
    mo = cand.iloc[0]
    motif = mo.motif
    used.add(motif)
    n_sites = sum(1 for _ in open(os.path.join(BASE, 'AMY_MC_footprints.bw', motif, 'beds', f'{motif}_all.bed')))
    ntf = int(ntf_map.get(motif, 1))

    fig, axes = plt.subplots(3, 1, figsize=(7.2, 7.4), sharex=True)
    for i, reg in enumerate(regions):
        mc, mw = read_signal(reg, motif)
        x = np.arange(len(mc)) - len(mc) // 2
        ax = axes[i]
        ax.plot(x, mc, color='#2c6fbb', lw=1.1, label='MC')
        ax.plot(x, mw, color='#c0392b', lw=1.1, label='MW', alpha=.9)
        ax.set_ylabel(f'{reg}\ncut signal', fontsize=8)
        ax.tick_params(labelsize=8)
        ax.spines[['top', 'right']].set_visible(False)
        if i == 0:
            ax.legend(fontsize=7.5, frameon=False, ncol=2, loc='upper right')
    axes[-1].set_xlabel('bp from motif center', fontsize=8)
    shared = f'; motif shared by {ntf} TFs' if ntf > 1 else ''
    fig.suptitle(f'{name}  —  {motif}\n{n_sites} motif sites; OR={mo.or_fisher:.2f}, FDR={mo.fdr_fisher:.1e}, mean bf={mo.mean_bf:.3f}; {int(r.n_enriched)}/{int(r.n_motifs)} motifs enriched{shared}',
                 fontsize=9.5)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    stem = f"{name.replace('/','_')}__{motif}"
    fig.savefig(os.path.join(FP, stem + '.pdf'), bbox_inches='tight')
    fig.savefig(os.path.join(FP, stem + '.png'), dpi=160, bbox_inches='tight')
    plt.close(fig)
    rows.append(dict(tf=name, n_enriched=int(r.n_enriched), n_motifs=int(r.n_motifs),
                     motif=motif, n_sites=n_sites, or_fisher=mo.or_fisher,
                     fdr_fisher=mo.fdr_fisher, mean_bf=mo.mean_bf, n_tf_shared=ntf,
                     pdf=stem + '.pdf'))
    print('plotted', name, motif)

pd.DataFrame(rows).to_csv(os.path.join(FP, 'index.tsv'), sep='\t', index=False)
print('done', len(rows))
