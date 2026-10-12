"""Parse immutable XML and record auditable, provisional screening flags."""
import csv, json, re, pathlib, collections, xml.etree.ElementTree as ET
ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW = ROOT / 'output/01_检索与研究方案/raw'
OUT = ROOT / 'output/02_PubMed初步分析'
OUT.mkdir(exist_ok=True)
def txt(e): return '' if e is None else ''.join(e.itertext()).strip()
def csvout(path, rows, fields=None):
    rows=list(rows)
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0])); w.writeheader(); w.writerows(rows)
records=[]
for source in sorted(RAW.glob('pubmed_*.xml')):
    for node in ET.parse(source).getroot().findall('PubmedArticle'):
        c=node.find('MedlineCitation'); a=c.find('Article'); journal=a.find('Journal')
        year=txt(journal.find('JournalIssue/PubDate/Year'))
        date_text=txt(journal.find('JournalIssue/PubDate/MedlineDate'))
        if not year:
            m=re.search(r'\b(19|20)\d{2}\b',date_text); year=m.group() if m else ''
        authors=[]; affiliations=[]; author_affiliations=[]
        for au in a.findall('AuthorList/Author'):
            name=(txt(au.find('LastName'))+', '+(txt(au.find('ForeName')) or txt(au.find('Initials')))).strip(', ')
            name=name or txt(au.find('CollectiveName'))
            if name: authors.append(name)
            aff=[txt(e) for e in au.findall('AffiliationInfo/Affiliation')]
            affiliations.extend(aff); author_affiliations.append({'author':name,'affiliations':aff})
        abstract='\n'.join((e.get('Label','')+': ' if e.get('Label') else '')+txt(e) for e in a.findall('Abstract/AbstractText'))
        keywords=sorted(set(txt(e) for e in c.findall('KeywordList/Keyword') if txt(e)))
        mesh=[txt(e) for e in c.findall('MeshHeadingList/MeshHeading/DescriptorName')]
        ids={e.get('IdType'):txt(e) for e in node.findall('PubmedData/ArticleIdList/ArticleId')}
        types=[txt(e) for e in a.findall('PublicationTypeList/PublicationType')]
        r={'pmid':txt(c.find('PMID')),'doi':ids.get('doi','').lower(),'pmc':ids.get('pmc',''),
           'title':txt(a.find('ArticleTitle')),'year':year,'pubdate_text':ET.tostring(journal.find('JournalIssue/PubDate'),encoding='unicode'),
           'journal':txt(journal.find('Title')),'issn':txt(journal.find('ISSN')),
           'volume':txt(journal.find('JournalIssue/Volume')),'issue':txt(journal.find('JournalIssue/Issue')),
           'pages':txt(a.find('Pagination/MedlinePgn')) or ('–'.join(filter(None,[txt(a.find('Pagination/StartPage')),txt(a.find('Pagination/EndPage'))]))) or txt(a.find("ELocationID[@EIdType='pii']")),
           'authors':authors,'affiliations':sorted(set(affiliations)), 'author_affiliations':author_affiliations,
           'keywords':keywords,'mesh':mesh,'abstract':abstract,'languages':[txt(e) for e in a.findall('Language')],
           'publication_types':types,'source_file':source.name,'url':'https://pubmed.ncbi.nlm.nih.gov/'+txt(c.find('PMID'))+'/',
           'references_available':len(node.findall('PubmedData/ReferenceList/Reference'))}
        r['missing_fields']=[k for k in ('year','authors','affiliations','keywords') if not r[k]]
        r['type_flags']=[p for p in types if p in ['Case Reports','Editorial','Letter','Comment','Published Erratum','News','Congress','Retracted Publication']]
        r['screen_status']='pending_relevance_review'
        r['primary_exclusion']=''
        # Mutually exclusive hierarchical reasons; counts of raw tags also retained separately.
        for p in ['Case Reports','Editorial','Letter','Comment','Published Erratum','News','Congress','Retracted Publication']:
            if p in types: r['primary_exclusion']=p; break
        if not r['primary_exclusion'] and 'eng' not in r['languages']: r['primary_exclusion']='Non-English'
        if r['primary_exclusion']: r['screen_status']='provisional_exclude'
        r['relevance_decision']='pending'
        records.append(r)
for source in sorted(RAW.glob('pubmed_*.xml')):
    for node in ET.parse(source).getroot().findall('PubmedBookArticle'):
        c=node.find('BookDocument'); pmid=txt(c.find('PMID'))
        records.append({'pmid':pmid,'doi':'','pmc':'','title':txt(c.find('ArticleTitle')) or txt(c.find('Book/BookTitle')),
            'year':txt(c.find('Book/PubDate/Year')),'pubdate_text':'','journal':txt(c.find('Book/BookTitle')),
            'issn':'','authors':[],'affiliations':[],'author_affiliations':[],'keywords':[],'mesh':[],
            'abstract':txt(c.find('Abstract')),'languages':[txt(e) for e in c.findall('Language')],
            'publication_types':['Book chapter'],'source_file':source.name,
            'url':'https://pubmed.ncbi.nlm.nih.gov/'+pmid+'/','references_available':0,
            'missing_fields':[], 'type_flags':['Book chapter'],'screen_status':'provisional_exclude',
            'primary_exclusion':'Book chapter','relevance_decision':'not_applicable'})
expected=json.loads((RAW/'pubmed_esearch.json').read_text())['esearchresult']['idlist']
assert len(records)==len(expected) and set(r['pmid'] for r in records)==set(expected)
reviews_path=ROOT/'script/relevance_decisions.json'
reviews=json.loads(reviews_path.read_text()) if reviews_path.exists() else {}
for r in records:
    r['protocol_flag']=bool(re.search(r'\bprotocol\b',r['title'],re.I))
    if r['pmid'] in {'4414100','38188810','30377465'}:
        reason={'4414100':'Letter (title verified)','38188810':'Retraction notice','30377465':'Case report (two patients; abstract verified)'}[r['pmid']]
        r['primary_exclusion']=reason;r['screen_status']='provisional_exclude'
        r['manual_type_note']='题名/摘要复核发现文献类型标记未充分反映实际类型。'
    if r['pmid'] in reviews:
        d=reviews[r['pmid']]; r['relevance_decision']=d['decision'];r['relevance_reason']=d['reason']
        if not r['primary_exclusion']:
            if d['decision']=='exclude':r['primary_exclusion']='Unrelated';r['screen_status']='provisional_exclude'
            elif d['decision']=='include':
                r['screen_status']='pending_metadata_completion' if r['missing_fields'] else ('provisional_complete' if r['abstract'] else 'pending_fulltext_check')
(OUT/'records.json').write_text(json.dumps(records,ensure_ascii=False,indent=2))
flat=[]
for r in records:
    x={k:json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else v for k,v in r.items()}
    x['reviewer_1']='';x['reviewer_2']='';x['final_decision']='';flat.append(x)
fields=list(dict.fromkeys(k for r in flat for k in r))
csvout(OUT/'01_逐篇筛选与书目.csv',flat,fields)
summary={'retrieved':len(records),'english':sum('eng' in r['languages'] for r in records),
    'raw_type_counts':dict(collections.Counter(p for r in records for p in r['type_flags'])),
    'primary_exclusion_counts':dict(collections.Counter(r['primary_exclusion'] for r in records if r['primary_exclusion'])),
    'status_counts':dict(collections.Counter(r['screen_status'] for r in records)),
    'missing_field_counts_before_merge':dict(collections.Counter(p for r in records for p in r['missing_fields'])),
    'duplicate_pmid':len(records)-len(set(r['pmid'] for r in records)),
    'editorial_pubmed_tag':sum('Editorial' in r['publication_types'] for r in records),
    'protocol_title_count':sum(r['protocol_flag'] for r in records),
    'screening_note':'AI-assisted partial screening. Not independently reviewed by two researchers. Missing metadata exclusions deferred until WoS merge.'}
(OUT/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
print(json.dumps(summary,ensure_ascii=False,indent=2))
