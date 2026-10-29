"""Independent checks, run AFTER the complete R analysis; internal, not for public supplement."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from scipy.stats import false_discovery_control
root=Path(__file__).resolve().parents[3]
b=Path(__file__).resolve().parent
D=pd.read_csv(b/'20261024_03_all_measured_gene_DE.csv')
X=pd.read_csv(b/'20261024_14_patient_group_expression.csv',index_col=0)
M=pd.read_csv(b/'20261024_02_sample_manifest_all.csv')
G=pd.read_csv(b/'20261024_13_patient_group_samples.csv')
E=pd.read_csv(b/'20261024_04_candidate_set_enrichment.csv')
assert len(M)==92
assert M.GSM.is_unique and D.gene.is_unique and X.index.is_unique
assert set(M.loc[M.included_in_paired_analysis,'tissue_group'])=={'PTC','Normal thyroid'}
assert G.groupby('patient_id').tissue_group.nunique().eq(2).all()
assert G.patient_id.nunique()==24 and len(G)==48
assert sorted(G.groupby('tissue_group').block_n.sum().tolist())==[27,28]
cols=list(X.columns)
# R write.csv(check.names=FALSE) preserves patient|tissue names.
ptc={c.split('|')[0]:c for c in cols if c.endswith('|PTC')}
normal={c.split('|')[0]:c for c in cols if c.endswith('|Normal thyroid')}
assert set(ptc)==set(normal) and len(ptc)==24
patients=sorted(ptc)
mean_fc=(X[[ptc[p] for p in patients]].to_numpy()-X[[normal[p] for p in patients]].to_numpy()).mean(axis=1)
R=D.set_index('gene').loc[X.index]
max_fc_error=float(np.max(np.abs(mean_fc-R.logFC.to_numpy())))
assert max_fc_error<1e-8,max_fc_error
assert np.isfinite(D[['logFC','P.Value','adj.P.Val']].to_numpy()).all()
assert D['P.Value'].between(0,1).all() and D['adj.P.Val'].between(0,1).all()
bh=false_discovery_control(D['P.Value'].to_numpy(),method='bh')
max_bh_error=float(np.max(np.abs(bh-D['adj.P.Val'].to_numpy())))
assert max_bh_error<1e-10,max_bh_error
sig=(D['adj.P.Val']<.05)&(D.logFC.abs()>=1)
assert np.array_equal(sig,D.DE_threshold)
measured=set(D.gene); de=set(D.loc[sig,'gene'])
disease=set(pd.read_csv(root/'v20241015/input/GeneCards-SearchResults_Thyroid_tumor.csv')['Gene Symbol'])
checks=[]
for _,r in E.iterrows():
 bg=measured if r.background=='all_measured' else measured&disease
 targets=set(pd.read_csv(root/f'v20241015/output/02_common_gene_{r["set"]}.csv').Com_Gene)&bg
 dx=de&bg
 a=len(targets&dx); bb=len(targets-dx); c=len(dx-targets); dd=len(bg-targets-dx)
 odds,p=fisher_exact([[a,bb],[c,dd]],alternative='greater')
 assert len(bg)==r.background_n and len(targets)==r.set_n and a==r.DE_in_set
 assert np.isclose(p,r.p_value,rtol=1e-8,atol=1e-300),(p,r.p_value)
 checks.append({'set':r['set'],'background':r.background,'a':a,'b':bb,'c':c,'d':dd,'p_python':p})
assert np.allclose(false_discovery_control(E.p_value.to_numpy()),E.FDR_6tests.to_numpy(),atol=1e-12)
report={'passed':True,'paired_patients':24,'observations':48,'measured_genes':len(D),'DE_genes':int(sig.sum()),'up_DE':int((sig&(D.logFC>0)).sum()),'down_DE':int((sig&(D.logFC<0)).sum()),'max_abs_logFC_error':max_fc_error,'max_abs_BH_error':max_bh_error,'fisher_checks':checks,'prespecified_genes':D[D.gene.isin(['CYP3A4','CACNA1H','POLB'])].to_dict(orient='records')}
(b/'20261024_29_independent_review.json').write_text(json.dumps(report,indent=2,ensure_ascii=False))
print(json.dumps(report,indent=2,ensure_ascii=False))
