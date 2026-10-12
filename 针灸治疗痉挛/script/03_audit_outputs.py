"""Search-quality diagnostics only: these are NOT final bibliometric results."""
import os, pathlib, json, csv, collections, hashlib, datetime
ROOT=pathlib.Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/acupuncture_spasm_matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
RAW=ROOT/'output/01_检索与研究方案/raw'; OUT=ROOT/'output/02_PubMed初步分析'
rs=json.loads((OUT/'records.json').read_text()); summary=json.loads((OUT/'summary.json').read_text())
base=json.loads((RAW/'pubmed_esearch.json').read_text())['esearchresult']
expanded=json.loads((RAW/'pubmed_sensitivity_esearch.json').read_text())['esearchresult']
added=sorted(set(expanded['idlist'])-set(base['idlist']),key=int)
removed=set(base['idlist'])-set(expanded['idlist'])
assert not removed, 'The sensitivity query must retain all baseline search results.'
def table(name,rows,fields):
    with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)
table('02_排除原因初筛统计.csv',[{'reason':k,'count':v,'status':'provisional; partial screening'} for k,v in summary['primary_exclusion_counts'].items()],['reason','count','status'])
years=range(min(int(r['year']) for r in rs if r['year']),2027)
allc=collections.Counter(int(r['year']) for r in rs if r['year'])
engc=collections.Counter(int(r['year']) for r in rs if r['year'] and 'eng' in r['languages'])
annual=[{'year':y,'all_retrieved':allc[y],'english_retrieved':engc[y],'partial_year':y==2026} for y in years]
table('03_原始检索年度分布_非最终发文量.csv',annual,['year','all_retrieved','english_retrieved','partial_year'])
table('04_敏感性检索新增PMID_未筛选.csv',[{'pmid':i,'url':f'https://pubmed.ncbi.nlm.nih.gov/{i}/'} for i in added],['pmid','url'])
articles=[r for r in rs if 'Book chapter' not in r['publication_types']]
table('05_元数据缺项审计.csv',[{'field':k,'missing':sum(not r.get(k) for r in articles),'denominator':len(articles)} for k in ['year','authors','affiliations','keywords','abstract']],['field','missing','denominator'])
fig,ax=plt.subplots(figsize=(11,4.6)); ax.bar(list(years),[allc[y] for y in years],color='#bed1dd',label='All retrieved records (n=395)')
ax.plot(list(years),[engc[y] for y in years],color='#136978',lw=1.8,label='English-language records (n=257)')
ax.axvspan(2025.5,2026.5,color='#dcab64',alpha=.4);ax.set(xlabel='Publication year in PubMed record',ylabel='Retrieved records',title='PubMed search audit — BEFORE final eligibility screening')
ax.spines[['top','right']].set_visible(False);ax.legend(frameon=False,loc='upper left');fig.text(.5,.01,'2026 is incomplete (cutoff: 10 October). Not a final publication-trend figure; no regression fitted.',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.05,1,1))
for ext in ['png','svg']:fig.savefig(OUT/f'06_原始检索年度审计.{ext}',dpi=220)
plt.close(fig)
labels=['Author keywords','Affiliations']; vals=[sum(not r['keywords'] for r in articles),sum(not r['affiliations'] for r in articles)]
fig,ax=plt.subplots(figsize=(8,3.7));ax.barh(labels,vals,color=['#176b7e','#aa7544'])
for i,v in enumerate(vals):ax.text(v+3,i,f'{v} / {len(articles)} journal records',va='center',fontsize=10)
ax.set(xlim=(0,285),xlabel='Records missing the field',title='Missing metadata before cross-database enrichment')
ax.invert_yaxis();ax.spines[['top','right']].set_visible(False);fig.text(.5,.01,'Overlapping flags, not exclusion counts. WoS may supply missing metadata.',ha='center',fontsize=9);fig.tight_layout(rect=(0,.05,1,1))
for ext in ['png','svg']:fig.savefig(OUT/f'07_元数据完整性审计.{ext}',dpi=220)
plt.close(fig)
data={'baseline':int(base['count']),'sensitivity':int(expanded['count']),'additional':len(added),
    'sentinel_42483336_in_baseline':'42483336' in base['idlist'], 'sentinel_42483336_in_sensitivity':'42483336' in expanded['idlist'],
    'query_translation':expanded['querytranslation'],
    'retrieved_at':datetime.datetime.fromtimestamp((RAW/'pubmed_sensitivity_esearch.json').stat().st_mtime,datetime.timezone.utc).isoformat(),
    'audit_generated_at':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(OUT/'sensitivity_comparison.json').write_text(json.dumps(data,ensure_ascii=False,indent=2))
(RAW/'checksums_sha256.json').write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in RAW.iterdir() if p.is_file() and p.name!='checksums_sha256.json'},indent=2))
report=f'''# PubMed 检索与初筛审计（非最终研究结果）

检索截止日期：2026-10-10（发表日期）；实际采集时间见 raw/manifest.json。原式检索 {len(rs)} 条，已下载全部 XML，PMID 集合与检索结果完全一致。XML 包含 394 篇期刊条目和 1 条 NICE 证据报告/书籍条目。不能把下载数量当作最终纳入数量。

## 当前计数

英文标记记录：{summary['english']}；非英文原始记录：{len(rs)-summary['english']}。采用文献类型优先的互斥原因后，单独计入“非英文”的记录为 131 条，另有非英文记录已归入更早的类型排除原因。

病例报告：PubMed 标记 49 条，另发现 1 篇双病例报告，共暂记 50 条。信函：PubMed 标记 1 条，另发现题名以 Letter 开头的 1 条，共暂记 2 条。社论：PubMed Editorial 标记为 0；Comment 为 1，单列，不擅自归为社论。其他：勘误 2，已撤稿原文 4，撤稿通知 1，NICE 报告/书籍条目 1。与研究主题无关：本轮已记录 10 条，逐条依据保存在筛选表。这不是完成全部相关性复核后的最终数量。

暂排除 202 条；尚未暂排除 193 条。其中 182 条仍待相关性判断，11 条已初步标记相关（7 条基本字段齐全，4 条仍缺字段）。字段齐全不等于最终纳入；全部仍需与 WoS 合并及最终复核。26 条原始记录题名含 protocol，未擅自将其当作已完成研究或全部删除。最终纳入数、WoS 数量、跨库重复数均未确定。

## 检索敏感性

在保留原式结果的前提下，单独试检补充 Muscle Spasticity MeSH、spastic*、spasm*、dry needl*、fire needl* 和 electro-acupuncture 拼写变体。得到 {data['sensitivity']} 条，比原式多 {data['additional']} 条；新增项仅留存 PMID，不混入主数据。该差额不能直接解释为真正漏检的合格文献，仍需筛选。

相近的2026年研究 PMID 42483336：原式是否命中 = {data['sentinel_42483336_in_baseline']}；扩展式是否命中 = {data['sentinel_42483336_in_sensitivity']}。这提示原式需要正式的查全率核查。

## 元数据与适用范围

394条期刊记录中，作者关键词缺失186条、机构缺失50条；NICE报告条目不纳入期刊字段完整性分母。MeSH与作者关键词分开保存，不以MeSH自动替代作者关键词。机构缺失统计仅指整条记录没有机构，并不代表其余记录每位作者均有机构。

补全之前不执行“缺机构/关键词”最终排除。所有原始摘要仅用于研究内部筛选与追溯。公开共享应遵守数据库与原出版物许可。题录包含PubMed链接、DOI、文献类型、语言、作者、机构、作者关键词、MeSH及原始文件来源。没有假称完成双人独立筛选。

## 图件状态

06、07号图为检索质量审计，不能替代SCI正文的正式发文趋势图或知识图谱。没有输出虚构的国家排名、中心性、Q/S、聚类标签或突现强度。正式图表的编号、算法与缺失输入见01文件夹中的研究方案。
'''
(OUT/'00_初筛审计说明.md').write_text(report,encoding='utf-8')
print(json.dumps(data,ensure_ascii=False,indent=2))
