import numpy as np, pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
B='/data2st2/junyi/output/atac1112/tobiasbam/BULK/BIND_DMR_enriched/'
m=pd.read_csv(B+'motif_enrichment_DMR_vs_genome_GCmatched.tsv',sep='\t')
en=(m.fdr_fisher<0.05)&(m.or_fisher>1)
e=m[en]; n=m[~en]

fig,axes=plt.subplots(1,2,figsize=(13.2,5.4))

ax=axes[0]
ax.scatter(n.or_fisher,n.mean_bf,s=9,c='#c8c8c8',alpha=.65,lw=0,label=f'not enriched (n={len(n)})')
ax.scatter(e.or_fisher,e.mean_bf,s=13,c='#d1495b',alpha=.8,lw=0,label=f'enriched, FDR<0.05 (n={len(e)})')
ax.axhline(.05,ls='--',lw=.8,c='#666'); ax.axhline(.10,ls=':',lw=.9,c='#666')
ax.text(1.02,.047,'bf = 5%',fontsize=7,color='#555',va='top')
ax.text(1.02,.096,'bf = 10%',fontsize=7,color='#555',va='top')
ax.set_xscale('log'); ax.set_xlim(.9,4.5)
ax.set_xlabel('motif enrichment OR (DMR vs GC-matched genome)')
ax.set_ylabel('mean bound fraction across 3 regions')
ann=[('metacluster_123.1','Irx2-6',(30,-6)),('tfdimers_MD00261','Pdx1',(22,10)),
     ('metacluster_165.3','Homez',(14,10)),('metacluster_67.2','E2f2/3',(18,10)),
     ('jaspar_MA1106.1','Hif1a',(-6,14)),('taipale_tf_pairs_E2F3_HES7_NNNNGCGCSNNNNNCACGTGNN_CAP','E2f3/Hes7',(10,-16))]
for mo,lab,off in ann:
    r=m[m.motif==mo]
    if len(r)==0: continue
    r=r.iloc[0]
    ax.annotate(lab,(r.or_fisher,r.mean_bf),fontsize=7.5,xytext=off,textcoords='offset points',color='#222',
                arrowprops=dict(arrowstyle='-',lw=.55,color='#777'))
ax.legend(fontsize=8,frameon=False,loc='lower right')
ax.set_title('Footprint support vs sequence enrichment',fontsize=10)

ax=axes[1]
bands=[('bf > 10%',.10,9),('5-10%',.05,3),('bf <= 5%',-1,1)]
def count(df,lo,hi):
    return ((df.mean_bf>lo)&(df.mean_bf<=hi)).sum()
vals=np.array([[count(e,.10,9),count(e,.05,.10),count(e,-1,.05)],
               [count(n,.10,9),count(n,.05,.10),count(n,-1,.05)]])
x=np.arange(2); w=.55; colors=['#2a6f97','#89b0c9','#e0e0e0']
bot=np.zeros(2)
for j,(lab,_,_) in enumerate(bands):
    ax.bar(x,vals[:,j],w,bottom=bot,color=colors[j],edgecolor='white',label=lab)
    for i in range(2):
        if vals[i,j]>0:
            ax.text(x[i],bot[i]+vals[i,j]/2,str(vals[i,j]),ha='center',va='center',fontsize=9,
                    color='white' if j!=2 else '#333',fontweight='bold')
    bot+=vals[:,j]
ax.set_xticks(x); ax.set_xticklabels([f'enriched\n(n={len(e)})',f'not enriched\n(n={len(n)})'])
ax.set_ylabel('# motifs'); ax.legend(fontsize=8,frameon=False,loc='upper right')
ax.set_title('Bound-fraction distribution',fontsize=10)
ax.spines[['top','right']].set_visible(False)

fig.tight_layout()
fig.savefig('/tmp/enrich/footprint_support.png',dpi=170,bbox_inches='tight')
print(vals)
print('saved')
