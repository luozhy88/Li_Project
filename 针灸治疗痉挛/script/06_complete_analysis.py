"""Reproducible PubMed-only completion; preserves the original audit artifacts."""
from pathlib import Path
import json, subprocess, hashlib, datetime, xml.etree.ElementTree as ET
ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/'output/01_检索与研究方案/raw'
OUT=ROOT/'output/06_完整分析'
RAW=OUT/'raw'; RAW.mkdir(parents=True,exist_ok=True)
base=json.loads((OLD/'pubmed_esearch.json').read_text())['esearchresult']['idlist']
expanded=json.loads((OLD/'pubmed_sensitivity_esearch.json').read_text())['esearchresult']['idlist']
added=[x for x in expanded if x not in base]
for start in range(0,len(added),100):
    dest=RAW/f'additional_{start:04d}.xml'
    if dest.exists(): continue
    temp=dest.with_suffix('.tmp')
    subprocess.run(['curl','-sS','-L','--fail','--retry','3','--max-time','90','-G',
      'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi','--data-urlencode','db=pubmed',
      '--data-urlencode','retmode=xml','--data-urlencode','id='+','.join(added[start:start+100]),'-o',str(temp)],check=True)
    root=ET.parse(temp).getroot()
    actual=[x.text for x in root.findall('PubmedArticle/MedlineCitation/PMID')]+[x.text for x in root.findall('PubmedBookArticle/BookDocument/PMID')]
    assert set(actual)==set(added[start:start+100]), (len(actual),len(added[start:start+100]))
    temp.rename(dest)
    print('Downloaded',dest.name,len(actual),flush=True)
manifest={'retrieved_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),
          'archived_search_date':'2026-10-11','publication_cutoff':'2026-10-10','baseline':len(base),'expanded':len(expanded),'added':len(added),
          'files':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in list(OLD.glob('*.xml'))+list(RAW.glob('*.xml'))}}
(OUT/'source_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))

# Reuse the original XML field parser without executing its output-writing section.
source=(ROOT/'script/02_parse_screen.py').read_text().split('expected=json.loads')[0]
source=source.replace("sorted(RAW.glob('pubmed_*.xml'))",'sources')
ns={'__file__':str(ROOT/'script/02_parse_screen.py'),'sources':sorted(OLD.glob('pubmed_*.xml'))+sorted(RAW.glob('additional_*.xml'))}
exec(compile(source,'original_xml_parser','exec'),ns)
records=ns['records']; byid={r['pmid']:r for r in records}
assert len(records)==len(byid)==604 and set(byid)==set(expanded)
for path in ns['sources']:
    for node in ET.parse(path).getroot().findall('PubmedArticle'):
        c=node.find('MedlineCitation'); r=byid[c.findtext('PMID')]
        r['keyword_sources']=[{'owner':k.get('Owner',''), 'terms':[ns['txt'](e) for e in k.findall('Keyword')]} for k in c.findall('KeywordList')]
        r['keywords']=sorted({t for k in r['keyword_sources'] if k['owner']=='NOTNLM' for t in k['terms'] if t})
        r['keyword_definition']='PubMed NOTNLM keywords (publisher-supplied; not independently verified as author keywords)'
        r['dates']=[{'type':d.get('DateType'), 'value':'-'.join(d.findtext(t,'') for t in ['Year','Month','Day'])} for d in c.findall('Article/ArticleDate')]
        r['publication_status']=node.findtext('PubmedData/PublicationStatus','')
        r['comments_corrections']=[{'type':e.get('RefType'),'pmid':e.findtext('PMID','')} for e in c.findall('CommentsCorrectionsList/CommentsCorrections')]
        r['references']=[{'citation':ns['txt'](e.find('Citation')),'pmid':e.findtext("ArticleIdList/ArticleId[@IdType='pubmed']",'')} for e in node.findall('PubmedData/ReferenceList/Reference')]
for r in records:
    r['retrieval_group']='original' if r['pmid'] in base else 'additional'
(OUT/'all_records.json').write_text(json.dumps(records,ensure_ascii=False,indent=2))
print('Parsed',len(records),'records')
