#!/usr/bin/env python
"""Annotate the 1126 v10nr_clust motifs with their TF(s) using the pySCENIC/cisTarget
motif annotation table ('motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl').

Mapping rules
-------------
* The MEME motif ids use single-underscore names ('swissregulon_hs_CUX1'); the tbl
  uses double underscores ('swissregulon__hs__CUX1'), so both sides are normalised.
* 'direct'  = rows with motif_similarity_qvalue == 0 (the motif itself is annotated
  for that gene in mouse).
* 'ortholog'= no direct row, best (smallest qvalue) annotation(s) instead.
"""
import pandas as pd, numpy as np

TBL = '/data2st1/junyi/scenic/mouse/motif/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl'
BASE = '/data2st2/junyi/output/atac1112/tobiasbam/BULK/BIND_DMR_enriched'
MEME = f'{BASE}/motifs/v10nr_clust.from_hdf5_DMRpycisfimoenriched.named.meme'
ENR = f'{BASE}/motif_enrichment_DMR_vs_flank.tsv'

cols = ['motif_id', 'motif_name', 'motif_description', 'source_name', 'source_version',
        'gene_name', 'qvalue', 'similar_motif_id', 'similar_motif_description',
        'orthologous_identity', 'orthologous_gene_name', 'orthologous_species', 'description']
ann = pd.read_csv(TBL, sep='\t', comment='#', header=None, names=cols, dtype=str, low_memory=False)
ann['qvalue'] = pd.to_numeric(ann['qvalue'], errors='coerce')
ann['norm_id'] = ann.motif_id.str.replace('__', '_', regex=False)

ids = sorted({l.split()[1] for l in open(MEME) if l.startswith('MOTIF ')})
print('motifs:', len(ids))

by_id = {k: v for k, v in ann.groupby('norm_id')}
rows = []
for m in ids:
    sub = by_id.get(m)
    if sub is None:
        rows.append((m, '', 0, 'none', '', np.nan, ''))
        continue
    direct = sub[sub.qvalue == 0]
    members = ';'.join(sorted(set(direct.motif_name))) if len(direct) else ''
    if len(direct):
        genes = sorted(set(direct.gene_name.dropna().str.strip()))
        lvl = 'direct'
        tf = ';'.join(genes)
        q = 0.0
    else:
        qmin = sub.qvalue.min()
        best = sub[sub.qvalue == qmin]
        genes = sorted(set(best.gene_name.dropna().str.strip()))
        lvl = 'ortholog'
        tf = ';'.join(genes)
        q = float(qmin) if pd.notna(qmin) else np.nan
    rows.append((m, tf, len(tf.split(';')) if tf else 0, lvl, members, q,
                 ';'.join(sorted(set(sub.source_name)))))

out = pd.DataFrame(rows, columns=['motif', 'tf', 'n_tf', 'annot_level', 'member_motifs',
                                  'min_qvalue', 'sources'])
out.to_csv(f'{BASE}/motif_tf_annotation.tsv', sep='\t', index=False)
print('annot level:', out.annot_level.value_counts().to_dict())
print('motifs without any annotation:', int((out.annot_level == 'none').sum()))
print(out[out.motif.isin(['metacluster_67.2', 'c2h2_zfs_M0369', 'hocomoco_ZBT14_HUMAN.H11MO.0.C',
                          'metacluster_55.3', 'taipale_cyt_meth_GMEB1_NKACGTNNNNACGTMN_FL_meth_repr',
                          'homer_NANTGCSGCA_Zfp57'])][['motif', 'tf', 'n_tf', 'annot_level']].to_string(index=False))

if __import__('os').path.exists(ENR):
    enr = pd.read_csv(ENR, sep='\t')
    merged = enr.merge(out, on='motif', how='left')
    merged.to_csv(f'{BASE}/motif_enrichment_DMR_vs_flank_with_TF.tsv', sep='\t', index=False)
    print('wrote merged table', merged.shape)
