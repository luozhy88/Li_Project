"""Save a dated PubMed search snapshot without replacing the frozen cohort."""
from pathlib import Path
import json,subprocess,hashlib,datetime
R=Path(__file__).resolve().parents[1];O=R/'output/13_联网复核';D=O/'raw';D.mkdir(parents=True,exist_ok=True)
source=R/'output/01_检索与研究方案/raw/pubmed_esearch.json'
old=json.loads(source.read_text())['esearchresult'];query=old['querytranslation']
log=[];response=None
for route in ['configured_proxy','direct']:
 p=D/('pubmed_original_'+route+'.json')
 cmd=['curl','-fsSL','--retry','1','--retry-all-errors','--retry-delay','1','--max-time','12']
 if route=='direct':cmd+=['--noproxy','*']
 cmd+=['-G','https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi','--data-urlencode','db=pubmed','--data-urlencode','retmode=json','--data-urlencode','retmax=9999','--data-urlencode','term='+query,'-o',str(p)]
 r=subprocess.run(cmd,capture_output=True,text=True)
 entry={'route':route,'exit_code':r.returncode,'error':r.stderr.strip()}
 if r.returncode==0:
  try:
   d=json.loads(p.read_text());s=d['esearchresult'];int(s['count']);entry['count']=int(s['count']);response=s;entry['sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
  except (ValueError,KeyError) as e:entry['parse_error']=str(e)
 log.append(entry)
 if response is not None:break
summary={'checked_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'query':query,'publication_cutoff':'2026-10-10','old_original_count':int(old['count']),'attempts':log,'status':'search_snapshot_saved' if response else 'search_not_retrieved','note':'Current retrospective search; does not replace the archived cohort or represent indexing as of the publication cutoff.'}
if response is not None:
 summary.update(current_count=int(response['count']),added_pmids=sorted(set(response['idlist'])-set(old['idlist'])),removed_pmids=sorted(set(old['idlist'])-set(response['idlist'])),actual_query_translation=response.get('querytranslation'),warnings=response.get('warninglist'))
(O/'pubmed_refresh_audit.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
print(json.dumps(summary,ensure_ascii=False,indent=2))
