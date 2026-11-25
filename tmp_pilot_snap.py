import time, sys, os
import numpy as np
import snapatac2 as snap

DS = '/data2st1/junyi/output/atac0627/doublet_filtered.h5ads/_dataset.h5ads'
OUT = '/data2st1/junyi/output/atac1112/l1_tracks_snapatac'
BL = '/data2st1/junyi/ref/mm10-blacklist.v2.bed'

t0 = time.time()
print('reading', DS, flush=True)
d = snap.read_dataset(DS)
print('loaded', d.shape, f'{time.time()-t0:.1f}s', flush=True)

reg = np.array(d.obs['Region'])
con = np.array(d.obs['Condition'])
l1 = np.array(d.obs['celltype.L1_ct'])

def export(labels, group, prefix, outdir=OUT):
    os.makedirs(outdir, exist_ok=True)
    t = time.time()
    n = int(np.sum(labels == group))
    print(f'EXPORT {prefix} group={group}: n={n}', flush=True)
    res = snap.ex.export_coverage(
        d, groupby=labels.tolist(), selections=[group],
        bin_size=10, normalization='CPM', blacklist=BL,
        max_frag_length=2000, counting_strategy='fragment',
        out_dir=outdir, prefix=prefix, suffix='.bw',
        tempdir='/dev/shm/snaptmp', n_jobs=32,
    )
    for k, v in res.items():
        print(f'  -> {k}: {v} ({os.path.getsize(v)/1e9:.3f} GB) in {time.time()-t:.1f}s', flush=True)
    return res

labels = np.array(['DROP'] * d.n_obs, dtype=object)
labels[(reg == 'AMY') & (con == 'MC')] = 'AMY_MC'
labels[(reg == 'AMY') & (con == 'MC') & (l1 == 'Glut')] = 'AMY_MC_Glut'
print('group counts', {g: int((labels == g).sum()) for g in ['AMY_MC', 'AMY_MC_Glut']}, flush=True)

export(labels, 'AMY_MC_Glut', 'pilot_AMY_MC_Glut_')
export(labels, 'AMY_MC', 'pilot_AMY_MC_all_')
print('DONE total', f'{time.time()-t0:.1f}s', flush=True)
