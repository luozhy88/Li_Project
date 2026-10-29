from pathlib import Path
import pandas as pd,numpy as np,json
from scipy.stats import hypergeom
root=Path(__file__).resolve().parents[3];out=Path(__file__).resolve().parent
sets={g:set(pd.read_csv(root/f'v20241015/output/02_common_gene_{g}.csv').Com_Gene) for g in ['ZBM','XKC','XYT']}
counts=[];bg=[];checks=[]
for g in sets:
 d=pd.read_csv(root/f'v20241015/output/01_TCM_CID_compound_Gene_{g}.csv')
 linked=d.entrez_gene_symbol.notna() & d.entrez_gene_symbol.fillna('').str.strip().ne('')
 disease=d.entrez_gene_symbol.fillna('').map(lambda x:bool(set(x.split(';'))&sets[g]))
 counts.append({'group':g,'all_CIDs':d.CID.nunique(),'gene_linked_CIDs':d[linked].CID.nunique(),'disease_linked_CIDs':d[disease].CID.nunique()})
 for k in ['go','kegg']:
  e=pd.read_csv(root/f'v20241015/output/03_{k}_enrichment_results_{g}.csv')
  for ont,z in (e.groupby('ONTOLOGY') if k=='go' else [('KEGG',e)]):
   ns=z.GeneRatio.str.split('/').str[-1].astype(int).unique();Ns=z.BgRatio.str.split('/').str[-1].astype(int).unique()
   assert len(ns)==len(Ns)==1
   bg.append(dict(group=g,analysis=ont,input_annotated_genes=int(ns[0]),background_annotated_genes=int(Ns[0]),reported_terms=len(z)))
  ratios=e.GeneRatio.str.split('/').apply(lambda x:tuple(map(int,x)));background=e.BgRatio.str.split('/').apply(lambda x:tuple(map(int,x)))
  recomputed=np.array([hypergeom.sf(a-1,N,M,n) for (a,n),(M,N) in zip(ratios,background)])
  assert np.allclose(recomputed,e.pvalue,rtol=1e-8,atol=1e-14)
  e['raw_p_recomputed']=recomputed
  e['fold_enrichment']=[(a/n)/(M/N) for (a,n),(M,N) in zip(ratios,background)]
  e.to_csv(out/f'20261024_{g}_{k}_verified.csv',index=False)
  checks.append(dict(group=g,analysis=k,rows=len(e),max_abs_raw_p_difference=float(max(abs(recomputed-e.pvalue))),BH_status='Archived BH values retained; full tested hypothesis set unavailable for independent BH reconstruction.'))
pd.DataFrame(counts).to_csv(out/'20261024_compound_flow.csv',index=False)
pd.DataFrame(bg).to_csv(out/'20261024_enrichment_backgrounds.csv',index=False)
pd.DataFrame([{'group':g,'disease_genes':len(sets[g]),'outside_XYT':len(sets[g]-sets['XYT'])} for g in sets]).to_csv(out/'20261024_gene_nesting.csv',index=False)
(out/'20261024_statistical_checks.json').write_text(json.dumps(checks,indent=2))
print(counts);print('Verified',sum(x['rows'] for x in checks),'exported raw P values.')
