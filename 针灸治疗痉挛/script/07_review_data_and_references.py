"""Independent, non-destructive audit of manuscript data and reference metadata."""
from pathlib import Path
import collections
import csv
import hashlib
import json
import shutil
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'output/06_Medicine审稿修订'
RAW = ROOT / 'output/01_检索与研究方案/raw'
EVID = OUT / '核查证据'
EVID.mkdir(parents=True, exist_ok=True)

def txt(node):
    return '' if node is None else ''.join(node.itertext()).strip()

def dump(name, value):
    (EVID / name).write_text(json.dumps(value, ensure_ascii=False, indent=2))

def csvout(path, rows):
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

hashes = json.loads((RAW / 'checksums_sha256.json').read_text())
hash_check = {name: hashlib.sha256((RAW / name).read_bytes()).hexdigest() == value for name, value in hashes.items()}
assert all(hash_check.values()), hash_check
articles, books = [], []
for file in sorted(RAW.glob('pubmed_*.xml')):
    tree = ET.parse(file)
    articles.extend(tree.findall('PubmedArticle'))
    books.extend(tree.findall('PubmedBookArticle'))
ids = [txt(n.find('MedlineCitation/PMID')) for n in articles] + [txt(n.find('BookDocument/PMID')) for n in books]
base = set(json.loads((RAW / 'pubmed_esearch.json').read_text())['esearchresult']['idlist'])
expanded = set(json.loads((RAW / 'pubmed_sensitivity_esearch.json').read_text())['esearchresult']['idlist'])
assert set(ids) == base and len(ids) == len(base)
records = json.loads((ROOT / 'output/02_PubMed初步分析/records.json').read_text())
with (ROOT / 'output/02_PubMed初步分析/01_逐篇筛选与书目.csv').open(encoding='utf-8-sig', newline='') as f:
    screening = list(csv.DictReader(f))
assert {r['pmid'] for r in records} == {r['pmid'] for r in screening} == base
statuses = dict(collections.Counter(r['screen_status'] for r in screening))
assert statuses == dict(collections.Counter(r['screen_status'] for r in records))
final_nonempty = sum(bool(r.get('final_decision', '').strip()) for r in screening)
owners = collections.Counter()
keyword_rows = []
missing_keywords = missing_affiliations = english = both = missing_author_keywords = 0
corrected_keywords = {}
for node in articles:
    c = node.find('MedlineCitation'); a = c.find('Article')
    langs = [txt(e) for e in a.findall('Language')]
    english += 'eng' in langs
    kws = [txt(e) for e in c.findall('KeywordList/Keyword') if txt(e)]
    author_kws = [txt(e) for e in c.findall("KeywordList[@Owner='NOTNLM']/Keyword") if txt(e)]
    other_kws = [txt(e) for group in c.findall('KeywordList') if group.get('Owner','NLM') != 'NOTNLM' for e in group.findall('Keyword') if txt(e)]
    missing_author_keywords += not author_kws
    corrected_keywords[txt(c.find('PMID'))] = {'author_keywords':sorted(set(author_kws)), 'other_KeywordList_terms':sorted(set(other_kws))}
    affs = [txt(e) for e in a.findall('AuthorList/Author/AffiliationInfo/Affiliation') if txt(e)]
    missing_keywords += not kws
    missing_affiliations += not affs
    both += not kws and not affs
    kwowners = sorted(set(e.get('Owner', '') for e in c.findall('KeywordList') if any(txt(k) for k in e.findall('Keyword'))))
    owners.update(kwowners)
    keyword_rows.append({'PMID':txt(c.find('PMID')), 'KeywordList_Owner':';'.join(kwowners), 'keyword_count':len(kws), 'author_keyword_count':len(author_kws), 'affiliation_missing':not bool(affs)})
english += sum('eng' in [txt(e) for e in n.findall('BookDocument/Language')] for n in books)
audit = {
    'hash_checks':hash_check, 'raw_records':len(ids), 'journal_records':len(articles), 'book_records':len(books),
    'unique_pmid':len(set(ids)), 'english_tag':english, 'missing_any_KeywordList':missing_keywords,
    'missing_author_keywords_NOTNLM':missing_author_keywords,
    'missing_any_affiliation':missing_affiliations, 'missing_both':both,
    'keyword_owner_record_counts':dict(owners), 'expanded_records':len(expanded),
    'base_not_in_expanded':len(base - expanded), 'expanded_new':len(expanded - base),
    'screen_status_counts':statuses, 'final_decision_nonempty':final_nonempty,
    'reviewer_1_nonempty':sum(bool(r['reviewer_1'].strip()) for r in screening),
    'reviewer_2_nonempty':sum(bool(r['reviewer_2'].strip()) for r in screening),
    'near_neighbor_42483336_in_base':'42483336' in base,
    'near_neighbor_42483336_in_expanded':'42483336' in expanded,
    'actual_search_time':json.loads((RAW / 'manifest.json').read_text())['retrieved_at'],
}
assert (len(ids), len(articles), english, missing_keywords, missing_affiliations, len(expanded)) == (395,394,257,186,50,604)
dump('01_独立数据复算.json', audit)
dump('05_作者关键词与其他词源分离.json',corrected_keywords)
csvout(EVID / '02_关键词字段来源逐篇核查.csv', keyword_rows)
rows = [
    ('原式检索记录',len(ids),395,'pubmed_*.xml + pubmed_esearch.json','检索候选数；非最终纳入数'),
    ('期刊记录',len(articles),395,'PubmedArticle节点','缺项统计分母'),
    ('英文标记',english,395,'Article/Language及BookDocument/Language','语言字段；不表示合格纳入'),
    ('缺少KeywordList关键词',missing_keywords,394,'MedlineCitation/KeywordList','并非严格作者关键词缺失数'),
    ('缺少作者关键词',missing_author_keywords,394,'KeywordList[@Owner=NOTNLM]','修订稿使用的作者关键词缺项数'),
    ('缺少机构字段',missing_affiliations,394,'Author/AffiliationInfo/Affiliation','没有任何机构字段；不等于机构信息完整率'),
    ('扩展检索记录',len(expanded),604,'pubmed_sensitivity_esearch.json','候选集'),
    ('扩展新增',len(expanded-base),395,'扩展PMID集合减原式集合','相对原式增加；不是合格漏检率'),
    ('暂标排除',statuses['provisional_exclude'],395,'逐篇筛选CSV','初筛状态'),
    ('未暂标排除',395-statuses['provisional_exclude'],395,'逐篇筛选CSV','不是已纳入数'),
    ('最终决定已填',final_nonempty,395,'逐篇筛选CSV/final_decision','空白不等于最终纳入0篇'),
]
csvout(OUT / '04_正文数据来源核对.csv', [{'指标':a,'复算数值':b,'分母':c,'百分比':round(100*b/c,2),'本地来源':d,'使用范围':e} for a,b,c,d,e in rows])

# Live bibliographic checks are saved separately, never merged into the study corpus.
if not (EVID / '03_引用文献实时PubMed.xml').exists():
    shutil.copyfile('/tmp/medicine_reference_pubmed.xml', EVID / '03_引用文献实时PubMed.xml')
fresh = {}
for n in ET.parse(EVID / '03_引用文献实时PubMed.xml').findall('PubmedArticle'):
    c=n.find('MedlineCitation'); a=c.find('Article'); j=a.find('Journal')
    pmid=txt(c.find('PMID'))
    doi=next((txt(e) for e in n.findall('PubmedData/ArticleIdList/ArticleId') if e.get('IdType')=='doi'), '')
    authors=[(txt(e.find('LastName'))+' '+txt(e.find('Initials'))).strip() for e in a.findall('AuthorList/Author')]
    fresh[pmid]={'pmid':pmid,'doi':doi,'title':txt(a.find('ArticleTitle')),'authors':authors,
        'journal':txt(j.find('ISOAbbreviation')), 'year':txt(j.find('JournalIssue/PubDate/Year')),
        'volume':txt(j.find('JournalIssue/Volume')),'issue':txt(j.find('JournalIssue/Issue')),
        'pages':txt(a.find('Pagination/MedlinePgn')) or '-'.join(filter(None,[txt(a.find('Pagination/StartPage')),txt(a.find('Pagination/EndPage'))])),
        'publication_types':[txt(e) for e in a.findall('PublicationTypeList/PublicationType')],
        'corrections':[{'type':e.get('RefType'),'pmid':txt(e.find('PMID'))} for e in c.findall('CommentsCorrectionsList/CommentsCorrections')],
        'article_dates':[{'type':e.get('DateType'),'date':'-'.join(txt(e.find(x)) for x in ['Year','Month','Day'])} for e in a.findall('ArticleDate')],
        'url':f'https://pubmed.ncbi.nlm.nih.gov/{pmid}/'}
dump('04_PubMed引用题录.json',fresh)
for doi in ['10.1002/asi.20317','10.1002/asi.21309','10.1186/s13643-020-01542-z','10.2196/57830']:
    filename=EVID / ('crossref_' + doi.replace('/','_') + '.json')
    if not filename.exists():
        req=urllib.request.Request('https://api.crossref.org/works/'+urllib.parse.quote(doi,safe=''),headers={'User-Agent':'ManuscriptReferenceAudit/1.0'})
        with urllib.request.urlopen(req,timeout=45) as response:
            payload=json.load(response)
        filename.write_text(json.dumps(payload,ensure_ascii=False,indent=2))
        time.sleep(0.4)
    msg=json.loads(filename.read_text())['message']
    assert msg['DOI'].lower()==doi.lower()
    print('Verified Crossref:',doi,msg['title'][0],flush=True)
print(json.dumps(audit,ensure_ascii=False,indent=2))
