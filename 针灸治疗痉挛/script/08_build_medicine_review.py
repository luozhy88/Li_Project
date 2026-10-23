"""Build the reviewed Chinese working manuscript and an independent Medicine-style table."""
import csv
import json
import re
from pathlib import Path
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'output/06_Medicine审稿修订'
EVID=OUT/'核查证据'
fresh=json.loads((EVID/'04_PubMed引用题录.json').read_text())
audit=json.loads((EVID/'01_独立数据复算.json').read_text())
refs={}

def expand_pages(pages):
    m=re.fullmatch(r'(\d+)-(\d+)',pages)
    if m and len(m[2])<len(m[1]):
        return m[1]+'-'+m[1][:-len(m[2])]+m[2]
    return pages

pmids={'park':'25192034','lim':'25628750','cai':'28455191','xue':'36062002','sanchez':'29986902','mu':'35245503','chen':'35837181','li':'42527851','liu':'42483336'}
for key,pmid in pmids.items():
    refs[key]=dict(fresh[pmid])
    refs[key]['pages']=expand_pages(refs[key]['pages'])
    refs[key]['status']='PubMed实时题录匹配；未见明确撤稿类型或关联'
    if key=='lim':
        refs[key]['status']='存在作者单位勘误PMID 28479927；非撤稿；核查PMC更正正文'
    if key=='li':
        refs[key]['status']='PubMed与Wiley匹配；在线先发表，卷期页留空'
for key,doi in [('qi','10.2196/57830'),('citespace','10.1002/asi.20317'),('cocitation','10.1002/asi.21309'),('prisma','10.1186/s13643-020-01542-z')]:
    msg=json.loads((EVID/('crossref_'+doi.replace('/','_')+'.json')).read_text())['message']
    authors=[]
    for a in msg['author']:
        if not a.get('family'):continue
        initials=''.join(p[0] for p in re.split(r'[\s.\-]+',a.get('given','')) if p)
        authors.append(a['family'].replace('‐','-')+' '+initials)
    year=str(msg.get('published-print',msg['published'])['date-parts'][0][0])
    refs[key]={'doi':doi,'pmid':'','title':msg['title'][0].replace('‐','-').rstrip('.'), 'authors':authors,
        'journal':'J Am Soc Inf Sci Technol' if key in ['citespace','cocitation'] else msg['short-container-title'][0],
        'year':year,'volume':msg.get('volume',''),'issue':msg.get('issue',''),
        'pages':msg.get('page',msg.get('article-number','')),'url':'https://doi.org/'+doi,
        'status':'Crossref登记题录匹配；结合出版方/作者机构页面核查；未见明确撤稿更新标记'}
refs['mesh']={'web':True,'title':'Muscle Spasticity. MeSH Descriptor Data 2026. D009128.','authors':['National Library of Medicine'],'url':'https://meshb.nlm.nih.gov/record/ui?ui=D009128','status':'NLM官方记录核实'}
refs['pubmed']={'web':True,'title':'PubMed User Guide.','authors':['National Library of Medicine'],'url':'https://pubmed.ncbi.nlm.nih.gov/help/','status':'NLM官方指南核实'}
refs['keyword']={'web':True,'title':'PubMed XML: KeywordList element.','authors':['National Library of Medicine'],'url':'https://dtd.nlm.nih.gov/ncbi/pubmed/doc/out/250101/el-KeywordList.html','status':'NLM官方字段说明核实；新增词源依据'}

def format_ref(r):
    authors=r['authors']
    author_text=', '.join(authors if len(authors)<=6 else authors[:3]+['et al'])+'. '
    if r.get('web'):
        return author_text+r['title']+' Accessed October 14, 2026. '+r['url']
    title=r['title'].rstrip('.')+'. '
    if r['volume']:
        tail=f"{r['journal']}. {r['year']};{r['volume']}"+(f"({r['issue']})" if r['issue'] else '')+f":{r['pages']}. doi:{r['doi']}"
    else:
        tail=f"{r['journal']}. Published online July 29, 2026. doi:{r['doi']}"
    return author_text+title+tail

body=(ROOT/'script/manuscript_medicine_reviewed.md').read_text()
pattern=r'\[([a-z]+(?:,[a-z]+)*)\]'
order=[]
for match in re.finditer(pattern,body):
    for key in match[1].split(','):
        assert key in refs,key
        if key not in order:order.append(key)
assert set(order)==set(refs)
number={k:i+1 for i,k in enumerate(order)}

def cite_numbers(keys):
    nums=sorted(number[k] for k in keys.split(','))
    groups=[];start=last=nums[0]
    for n in nums[1:]+[None]:
        if n is not None and n==last+1:
            last=n;continue
        if last-start>=2:groups.append(f'{start}-{last}')
        else:groups.extend(str(k) for k in range(start,last+1))
        start=last=n
    return ','.join(groups)

def doc_base(header):
    doc=Document();sec=doc.sections[0]
    sec.page_width=Inches(8.5);sec.page_height=Inches(11)
    sec.top_margin=sec.bottom_margin=sec.left_margin=sec.right_margin=Inches(0.85)
    for name in ['Normal','Title','Heading 1','Heading 2','Caption','Header','Footer']:
        st=doc.styles[name];st.font.name='Times New Roman';st.font.size=Pt(11)
        st.font.color.rgb=RGBColor.from_string('000000')
        st._element.get_or_add_rPr().get_or_add_rFonts().set(qn('w:eastAsia'),'Noto Serif CJK SC')
        st.paragraph_format.line_spacing=1.5
        st.paragraph_format.space_after=Pt(5)
        st.paragraph_format.widow_control=True
    for name,size in [('Title',17),('Heading 1',13),('Heading 2',11.5)]:
        doc.styles[name].font.size=Pt(size);doc.styles[name].font.bold=True
        doc.styles[name].paragraph_format.keep_with_next=True
    sec.header.paragraphs[0].text=header
    doc.styles['Header'].font.size=Pt(8)
    doc.styles['Footer'].font.size=Pt(8)
    p=sec.footer.paragraphs[0];p.alignment=2
    f=OxmlElement('w:fldSimple');f.set(qn('w:instr'),'PAGE');p._p.append(f)
    return doc

def add_cited(p,text):
    cursor=0
    for match in re.finditer(pattern,text):
        p.add_run(text[cursor:match.start()])
        run=p.add_run(cite_numbers(match[1]));run.font.superscript=True
        cursor=match.end()
    p.add_run(text[cursor:])

doc=doc_base('针灸相关干预与痉挛｜Medicine审稿修订版｜正式分析待完成')
for chunk in body.strip().split('\n\n'):
    chunk=chunk.strip()
    if chunk.startswith('# '):doc.add_paragraph(chunk[2:],'Title')
    elif chunk.startswith('## '):
        p=doc.add_paragraph(chunk[3:],'Heading 1')
        if chunk.startswith('## 摘要'):p.paragraph_format.page_break_before=True
    elif chunk.startswith('### '):doc.add_paragraph(chunk[4:],'Heading 2')
    else:add_cited(doc.add_paragraph(),chunk.replace('\n',' '))
p=doc.add_paragraph('参考文献','Heading 1');p.paragraph_format.page_break_before=True
for key in order:
    r=refs[key];text=f'{number[key]}. '+format_ref(r)
    p=doc.add_paragraph();p.paragraph_format.left_indent=Inches(.25);p.paragraph_format.first_line_indent=Inches(-.25)
    p.paragraph_format.keep_together=True
    p.paragraph_format.space_after=Pt(0)
    if not r.get('web'):
        before,after=text.split(r['journal'],1)
        p.add_run(before);p.add_run(r['journal']).italic=True;p.add_run(after)
    else:p.add_run(text)
    for run in p.runs:run.font.size=Pt(10)
doc.core_properties.title=body.splitlines()[0][2:]
doc.core_properties.subject='Medicine格式中文工作稿；核对本地数据及真实引文；非完整投稿终稿'
dest=OUT/'01_针灸治疗痉挛_Medicine审稿修订稿.docx'
doc.save(dest)
md=re.sub(pattern,lambda m:'['+cite_numbers(m[1])+']',body)
(OUT/'01_针灸治疗痉挛_Medicine审稿修订稿.md').write_text(md+'\n\n## 参考文献\n\n'+'\n\n'.join(f'{number[k]}. '+format_ref(refs[k]) for k in order)+'\n')

# Table is editable and separate, in line with the journal's submission instructions.
table_doc=doc_base('表1｜当前检索记录检查｜不是最终纳入样本统计')
table_doc.add_paragraph('表1 当前检索记录及字段检查结果','Heading 1')
table=table_doc.add_table(rows=1,cols=4);table.style='Table Grid'
for cell,value in zip(table.rows[0].cells,['指标','记录数','分母','比例（%）']):cell.text=value
table_rows=[('原式检索记录',395,395),('其中期刊记录',394,395),('具有英文语言标记',257,395),('具有作者关键词',207,394),('缺少作者关键词',187,394),('无任何关键词字段',186,394),('无机构字段',50,394),('暂标排除',202,395),('尚未暂标排除',193,395),('扩展检索候选记录',604,None),('扩展新增候选记录',209,None)]
for title,count,denom in table_rows:
    for cell,value in zip(table.add_row().cells,[title,str(count),str(denom) if denom else '不适用',f'{count/denom*100:.2f}' if denom else '不适用']):cell.text=value
for row in table.rows:
    for cell in row.cells:
        for p in cell.paragraphs:
            p.paragraph_format.space_after=Pt(3)
            for run in p.runs:run.font.size=Pt(10)
trpr=table.rows[0]._tr.get_or_add_trPr();repeat=OxmlElement('w:tblHeader');trpr.append(repeat)
table_doc.add_paragraph('注：数据来自本项目2026年10月11日保存的检索记录及本轮独立复算，非最终纳入集。作者关键词按PubMed XML中NOTNLM来源识别。1条记录仅含旧版NLM索引词，故187条缺少作者关键词与186条无任何关键词字段是不同指标；机构与关键词缺项可重叠。604条扩展记录包含全部395条原式记录，209条为集合差集。未对扩展候选作完整资格筛选。')
table_doc.save(OUT/'05_表1_本项目检索与字段核查.docx')

oldnum={'mesh':1,'park':2,'lim':3,'cai':4,'xue':5,'sanchez':6,'mu':7,'chen':8,'li':9,'liu':10,'citespace':11,'cocitation':12,'prisma':13,'pubmed':14,'qi':16,'keyword':'新增'}
audits=[]
for k in order:
    r=refs[k]
    audits.append({'原编号':oldnum[k],'修订编号':number[k],'题名':r['title'],'DOI':r.get('doi',''),'PMID':r.get('pmid',''),'来源链接':r['url'],'核查结论':r['status'],'处理':'保留并按Medicine格式重排' if k!='keyword' else '新增官方词源定义','核查日期':'2026-10-14'})
audits.append({'原编号':15,'修订编号':'移出正文','题名':'Web of Science: Export Records','DOI':'','PMID':'','来源链接':'https://webofscience.help.clarivate.com/en-us/Content/export-records.htm','核查结论':'旧链接跳转到官方帮助中心；本轮无法从原URL直接核对原导出说明','处理':'从正式引文移除；保留历史资源记录','核查日期':'2026-10-14'})
audits.append({'原编号':'原3的关联更正','修订编号':'仅审计','题名':fresh['28479927']['title'],'DOI':fresh['28479927']['doi'],'PMID':'28479927','来源链接':'https://pmc.ncbi.nlm.nih.gov/articles/PMC5396446/','核查结论':'第五作者单位更正；非撤稿、非疗效数据更正','处理':'在核查报告说明；不为凑数增加正文引文','核查日期':'2026-10-14'})
with (OUT/'03_逐条参考文献真实性核查.csv').open('w',encoding='utf-8-sig',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=list(audits[0]));writer.writeheader();writer.writerows(audits)
(OUT/'03_参考文献核查索引.md').write_text('# 参考文献核查索引\n\n题录核查不等于全文质量评价。原始响应见“核查证据”目录。\n\n'+'\n\n'.join(f"修订[{number[k]}]（原{oldnum[k]}）：{format_ref(refs[k])}\n\n核查：{refs[k]['status']}。来源：{refs[k]['url']}" for k in order)+'\n\n原[15]旧版Clarivate链接跳转，已移出正式引文。Lim等原[3]关联更正见PMID 28479927，只更正作者单位。\n')

# Full search details, copied from the frozen strategy rather than invented at revision time.
plan=(ROOT/'output/01_检索与研究方案/01_研究方案与图表执行规范.md').read_text()
search=plan.split('## 二、检索与数据来源')[1].split('### 2.3 需要的WoS导出')[0]
supp='# 补充材料S1：已保存检索策略\n\n本材料中的395/604条结果来自原保存时点，未在本轮重新检索替换。主检索策略尚需在双库正式分析前冻结。\n'+search
(OUT/'06_补充材料S1_检索策略.md').write_text(supp)
sd=doc_base('补充材料S1｜已保存检索策略')
in_code=False
for line in supp.splitlines():
    if line.startswith('```'):in_code=not in_code;continue
    if not line.strip():continue
    if line.startswith('# '):sd.add_paragraph(line[2:],'Title')
    elif line.startswith('### '):sd.add_paragraph(line[4:],'Heading 2')
    else:
        p=sd.add_paragraph(line)
        if in_code:
            p.paragraph_format.line_spacing=1.0
            for run in p.runs:run.font.size=Pt(10)
sd.save(OUT/'06_补充材料S1_检索策略.docx')

checked=Document(dest)
alltext='\n'.join(p.text for p in checked.paragraphs)
assert len(refs)==16 and len(order)==16
assert checked.sections[0].page_width==Inches(8.5)
assert checked.styles['Normal'].paragraph_format.line_spacing==1.5
assert len(checked.inline_shapes)==0
assert '187条（47.46%）' in alltext
assert '曹老师' not in alltext and '本节对应' not in alltext
assert not re.search(pattern,alltext)
assert audit['missing_author_keywords_NOTNLM']==187
qa={'reference_count':16,'journal_articles':13,'official_web_resources':3,'first_appearance_sequence':number,
    'page_size':'US Letter','line_spacing':1.5,'body_font':'Times New Roman / Noto Serif CJK SC, 11 pt',
    'images_embedded':0,'table_separate':True,'superscript_citation_runs':sum(bool(r.font.superscript) for p in checked.paragraphs for r in p.runs),
    'body_chinese_characters':len(re.findall(r'[\u4e00-\u9fff]',body)), 'formal_analysis_completed':False,
    'author_keyword_missing_corrected':187,'original_data_hashes_unchanged':all(audit['hash_checks'].values())}
(OUT/'核查证据/06_Word结构与引文检查.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2))
print(dest);print(json.dumps(qa,ensure_ascii=False))
