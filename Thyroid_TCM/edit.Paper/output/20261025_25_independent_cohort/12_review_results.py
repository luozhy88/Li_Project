from pathlib import Path
import pandas as pd,numpy as np,json,gzip,io
from scipy.stats import fisher_exact,spearmanr
p=Path(__file__).parent
meta=pd.read_csv(p/'11_verified_paired_manifest.csv',dtype={'patient':str})
expr=pd.read_csv(p/'GSE33630_14_patient_group_expression.csv',index_col=0)
de=pd.read_csv(p/'GSE33630_03_all_measured_gene_DE.csv',index_col=0)
ids=sorted(meta.patient.unique());diff=np.column_stack([expr[x+'|PTC']-expr[x+'|Normal thyroid'] for x in ids]);manual=diff.mean(axis=1)
assert len(ids)==44 and diff.shape==(20848,44)
err=float(np.max(abs(manual-de.loc[expr.index,'logFC'].values)));assert err<1e-9
pvals=de['P.Value'].values;order=np.argsort(pvals);bh=np.empty(len(pvals));bh[order]=np.minimum(1,np.minimum.accumulate((pvals[order]*len(pvals)/np.arange(1,len(pvals)+1))[::-1])[::-1]);assert np.max(abs(bh-de['adj.P.Val']))<1e-10
sets={k:set(pd.read_csv('v20241015/output/02_common_gene_'+k+'.csv').Com_Gene) for k in ['XYT','XKC','ZBM']}
dis=set(pd.read_csv('v20241015/input/GeneCards-SearchResults_Thyroid_tumor.csv')['Gene Symbol']);measured=set(de.index);deg=set(de.index[(de['adj.P.Val']<.05)&(de.logFC.abs()>=1)])
e=pd.read_csv(p/'GSE33630_enrichment_with_CI.csv');pv=[]
for _,r in e.iterrows():
 bg=measured if r.background=='all_measured' else measured&dis;s=sets[r['set']]&bg;d=deg&bg;a=len(s&d);b=len(s)-a;c=len(d)-a;dd=len(bg)-a-b-c
 v=fisher_exact([[a,b],[c,dd]],alternative='greater').pvalue;assert np.isclose(v,r.p_value,rtol=1e-8,atol=1e-100);assert a==r.DE_in_set;pv.append(v)
pv=np.array(pv);order=np.argsort(pv);q=np.empty(6);q[order]=np.minimum(1,np.minimum.accumulate((pv[order]*6/np.arange(1,7))[::-1])[::-1]);assert np.allclose(q,e.FDR_6tests,rtol=1e-8,atol=1e-100)
summary={'status':'PASS','paired_logFC_max_abs_error':err,'measured_genes':len(measured),'DE_genes':len(deg),'up':int(((de.logFC>=1)&(de['adj.P.Val']<.05)).sum()),'down':int(((de.logFC<=-1)&(de['adj.P.Val']<.05)).sum()),'sample_observations':expr.shape[1],'paired_patients':len(ids),'BH_and_six_Fisher_tests_independently_recomputed':True}
(p/'13_independent_numerical_review.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
