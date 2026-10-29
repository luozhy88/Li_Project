from pathlib import Path
import pandas as pd,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path(__file__).parent
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False})
d=pd.read_csv(p/'GSE33630_03_all_measured_gene_DE.csv');e=pd.read_csv(p/'GSE33630_enrichment_with_CI.csv');c=pd.read_csv(p/'GSE33630_cross_cohort_all_genes.csv')
fig,ax=plt.subplots(1,3,figsize=(13.5,4.3),gridspec_kw={'width_ratios':[1,1.15,1]})
mask=d.DE_threshold
ax[0].scatter(d.logFC,-np.log10(d['adj.P.Val'].clip(lower=1e-100)),s=3,c=np.where(mask,'#ac3939','#bbbbbb'),alpha=.6,rasterized=True)
ax[0].axvline(-1,color='gray',ls='--',lw=.7);ax[0].axvline(1,color='gray',ls='--',lw=.7);ax[0].axhline(-np.log10(.05),color='gray',ls='--',lw=.7)
ax[0].set(xlabel='log2 fold change (PTC vs normal)',ylabel='−log10 BH-adjusted P',title='A  GSE33630: 44 matched pairs')
ax[0].text(.04,.94,'745 DE genes\n407 up / 338 down',transform=ax[0].transAxes,va='top')
colors={'XYT':'#7851a9','XKC':'#208b8b','ZBM':'#dc8a29'}
for i,(_,r) in enumerate(e.iterrows()):
 y=5-i;ax[1].errorbar(r.odds_ratio,y,xerr=[[r.odds_ratio-r.CI_low],[r.CI_high-r.odds_ratio]],fmt='o',color=colors[r['set']],capsize=3,markersize=5)
 ax[1].text(4,y,f'q={r.FDR_6tests:.2g}',va='center',fontsize=8)
ax[1].set_yticks(range(6),[f"{r['set']} / {'all' if r.background=='all_measured' else 'disease'}" for _,r in e.iloc[::-1].iterrows()]);ax[1].axvline(1,ls='--',color='gray',lw=.8);ax[1].set_xscale('log');ax[1].set_xlim(.05,12);ax[1].set_ylim(-.8,5.8);ax[1].set(xlabel='Odds ratio (exact 95% CI)',title='B  Candidate-set enrichment')
ax[2].scatter(c.logFC_first,c.logFC_second,s=3,c='#a5a5a5',alpha=.3,rasterized=True)
f=c[c.DE_threshold_first];ax[2].scatter(f.logFC_first,f.logFC_second,s=5,c='#2878a5',alpha=.65,rasterized=True)
ax[2].axhline(0,color='gray',lw=.7);ax[2].axvline(0,color='gray',lw=.7);ax[2].set(xlabel='GSE60542 log2 fold change',ylabel='GSE33630 log2 fold change',title='C  Cross-cohort agreement')
ax[2].text(.04,.97,'Spearman ρ = 0.786\n664/666 directions agree\n536/666 replicate the DE threshold',transform=ax[2].transAxes,va='top',fontsize=8)
fig.tight_layout();fig.savefig(p/'16_Supplementary_Figure_S3.png',dpi=300,bbox_inches='tight');fig.savefig(p/'16_Supplementary_Figure_S3.pdf',bbox_inches='tight');plt.close(fig)
