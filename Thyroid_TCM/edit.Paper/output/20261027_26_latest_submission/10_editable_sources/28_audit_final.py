"""Read-only scientific checks; writes audit artifacts, never analysis inputs."""
from pathlib import Path
import json, re, subprocess
from urllib.parse import urlencode
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr, false_discovery_control
from scipy.stats.contingency import odds_ratio

pkg = Path(__file__).resolve().parents[1]
root = pkg.parent
audit = pkg / '11_final_review'
audit.mkdir(exist_ok=True)
md = (pkg/'01_Manuscript_EN.md').read_text()
ids = re.findall(r'\[PMID: (\d+)\]', md)
if '25417114' not in ids:
    ids.append('25417114')
url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?' + urlencode({
    'query': '('+' OR '.join('EXT_ID:'+i for i in ids)+') AND SRC:MED',
    'format':'json', 'resultType':'core', 'pageSize':100})
cache = audit/'reference_metadata_europepmc.json'
raw = cache.read_text() if cache.exists() and cache.stat().st_size else subprocess.check_output(['curl','-fsSL','--retry','2','--max-time','50',url], text=True)
(audit/'reference_metadata_europepmc.json').write_text(raw)
records = {r['id']:r for r in json.loads(raw)['resultList']['result']}
assert set(ids) == set(records), (set(ids)-set(records))
rows=[]
for i, pmid in enumerate(ids,1):
    r=records[pmid]; j=r['journalInfo']
    rows.append(dict(number=i, PMID=pmid, title=r['title'], authors=r['authorString'],
                     journal=j['journal']['medlineAbbreviation'],year=r['pubYear'],
                     volume=j.get('volume'),issue=j.get('issue'),pages=r.get('pageInfo'),
                     DOI=r.get('doi'),pubmed_url='https://pubmed.ncbi.nlm.nih.gov/'+pmid+'/',
                     publication_types='; '.join(r.get('pubTypeList',{}).get('pubType',[])),
                     correction_links=json.dumps(r.get('commentCorrectionList',{}))))
pd.DataFrame(rows).to_csv(audit/'reference_verification.csv',index=False)
print('REFERENCE METADATA',len(rows))

x = pd.ExcelFile(pkg/'05_Supplementary_Tables.xlsx')
sets={s:set(pd.read_excel(x,sheet_name=sheet)['Com_Gene']) for s,sheet in
      [('ZBM','T07_overlap_genes_ZBM'),('XKC','T08_overlap_genes_XKC'),('XYT','T09_overlap_genes_XYT')]}
disease=set(pd.read_csv(root.parent.parent/'v20241015/input/GeneCards-SearchResults_Thyroid_tumor.csv')['Gene Symbol'])
checks={}
for cohort,sheet,esheet in [('GSE60542','T14_GSE60542_all_gene_DE','T20_GSE60542_enrichment_CI'),('GSE33630','T23_GSE33630_all_gene_DE','T24_GSE33630_set_enrichment')]:
    d=pd.read_excel(x,sheet_name=sheet); e=pd.read_excel(x,sheet_name=esheet)
    e=e.rename(columns={'ci_low':'CI_low','ci_high':'CI_high','CI95_lower':'CI_low','CI95_upper':'CI_high','BH_adjusted_P_six_tests':'FDR_6tests'})
    measured=set(d.gene); de=set(d.loc[(d['adj.P.Val']<.05)&(d.logFC.abs()>=1),'gene'])
    np.testing.assert_allclose(false_discovery_control(d['P.Value']),d['adj.P.Val'],rtol=1e-10,atol=1e-14)
    ps=[]
    for _,r in e.iterrows():
        bg=measured if r.background=='all_measured' else measured & disease
        sg=sets[r['set']] & bg; dg=de & bg
        a=len(sg&dg); b=len(sg)-a; c=len(dg)-a; dd=len(bg)-a-b-c
        table=[[a,b],[c,dd]]
        ps.append(fisher_exact(table,alternative='greater').pvalue)
        o=odds_ratio(table,kind='conditional'); ci=o.confidence_interval()
        np.testing.assert_allclose(o.statistic,r.odds_ratio,rtol=1e-3)
        np.testing.assert_allclose([ci.low,ci.high],[r.CI_low,r.CI_high],rtol=1e-3)
    np.testing.assert_allclose(false_discovery_control(ps),e.FDR_6tests,rtol=1e-8)
    checks[cohort]={'genes':len(d),'DE':len(de),'up':int(((d.gene.isin(de))&(d.logFC>0)).sum()),'BH_Fisher_OR_CI':'PASS'}
cc=pd.read_excel(x,sheet_name='T25_cross_cohort_genes')
rep=cc.DE_threshold_first & cc.DE_threshold_second & cc.same_direction
checks['cross_cohort']={'rho':spearmanr(cc.logFC_first,cc.logFC_second).statistic,'same_direction_first_DE':int((cc.DE_threshold_first&cc.same_direction).sum()),'replicated':int(rep.sum()),'XYT_replicated':int((rep&cc.gene.isin(sets['XYT'])).sum()),'XKC_replicated':int((rep&cc.gene.isin(sets['XKC'])).sum())}
checks['set_relations']={s+'_subset_XYT':sets[s]<=sets['XYT'] for s in ['XKC','ZBM']}
checks['set_relations']['ZBM_subset_XKC']=sets['ZBM']<=sets['XKC']
t39=pd.read_excel(x,sheet_name='T39_replicated_candidates')
checks['T39_zero_adjusted_P_cells']=int((t39.filter(like='adj_P')==0).sum().sum())
checks['T39_top10']=t39.head(10).to_dict('records')
inventory=pd.read_csv(root.parent.parent/'v20241015/output/01_TCM_CID_compound_Gene_XYT.csv')
def gene_union(frame):
    return {g for v in frame.entrez_gene_symbol.dropna() for g in str(v).split(';') if g}
full=gene_union(inventory)
restricted=gene_union(inventory[inventory.herb_cn_name!='白毛夏枯草'])
assert full & disease == restricted & disease == sets['XYT']
checks['exclude_Bai_Mao_Xia_Ku_Cao']={'original_genes':len(full),'remaining_genes':len(restricted),'removed':sorted(full-restricted),'candidate_set_unchanged':True}
for cohort, folder, prefix, sheet in [
    ('GSE60542','20261024_25_transcriptome_validation','20261024','T14_GSE60542_all_gene_DE'),
    ('GSE33630','20261025_25_independent_cohort','GSE33630','T23_GSE33630_all_gene_DE')]:
    expr=pd.read_csv(root/folder/(prefix+'_14_patient_group_expression.csv'),index_col=0)
    normal=expr.filter(regex=r'\|Normal thyroid$')
    tumour=expr.filter(regex=r'\|PTC$')
    assert len(normal.columns)==len(tumour.columns) and len(normal.columns)>0
    delta=tumour.mean(axis=1)-normal.mean(axis=1)
    de=pd.read_excel(x,sheet_name=sheet).set_index('gene')
    err=float(np.max(np.abs(delta.loc[de.index]-de.logFC)))
    assert err < 1e-10
    checks[cohort]['paired_logFC_max_error']=err
strata=pd.read_excel(x,sheet_name='T27_exposure_sensitivity')
np.testing.assert_allclose(false_discovery_control(strata.p_value),strata.FDR_12tests,rtol=1e-8)
checks['exposure_sensitivity_BH']='PASS'
(audit/'numerical_audit.json').write_text(json.dumps(checks,indent=2))
print(json.dumps({k:v for k,v in checks.items() if k!='T39_top10'},indent=2))
