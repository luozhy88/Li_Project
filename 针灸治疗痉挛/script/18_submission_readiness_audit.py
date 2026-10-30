"""Audit archived data against the user's revised eligibility; never label this a dual-source final screen."""
from pathlib import Path
from collections import Counter
import csv,json,hashlib,datetime
R=Path(__file__).resolve().parents[1];O=R/'output/12_SCI投稿前核验';O.mkdir(exist_ok=True)
p=R/'output/09_审稿后重分析/screened_records.json';rs=json.loads(p.read_text())
required=['year','authors','affiliations','keywords']
def missing(r):return [k for k in required if not r.get(k)]
def write(name,rows,fields):
 with (O/name).open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
fields=['pmid','doi','title','retrieval_group','legacy_decision','legacy_exclusion','missing_core_fields','field_complete','has_abstract','has_reference_list','publication_types','year','affiliations','keywords','qualification_note']
write('01_604条记录_新标准字段核验.csv',[
 dict(pmid=r['pmid'],doi=r['doi'],title=r['title'],retrieval_group=r['retrieval_group'],legacy_decision=r['decision'],legacy_exclusion=r['exclusion_reason'],missing_core_fields='; '.join(missing(r)),field_complete=not missing(r),has_abstract=bool(r['abstract']),has_reference_list=bool(r.get('references')),publication_types='; '.join(r['publication_types']),year=r['year'],affiliations=' | '.join(r['affiliations']),keywords=' | '.join(r['keywords']),qualification_note='仅核验旧档案字段；未完成WoS补全、日期复核和新标准最终筛选') for r in rs],fields)
eligible=[r for r in rs if r['decision']=='include'];complete=[r for r in eligible if not missing(r)]
write('02_旧纳入但核心字段缺失_待跨库补全.csv',[{'pmid':r['pmid'],'doi':r['doi'],'title':r['title'],'missing_fields':'; '.join(missing(r)),'retrieval_group':r['retrieval_group'],'year':r['year']} for r in eligible if missing(r)],['pmid','doi','title','missing_fields','retrieval_group','year'])
reasons=Counter(r['exclusion_reason'] for r in rs if r['decision']=='exclude')
write('03_旧筛选互斥主原因.csv',[{'reason':k,'count':v,'scope':'旧604条PubMed扩展队列，非双库最终结果'} for k,v in sorted(reasons.items())],['reason','count','scope'])
# Separate all type labels from mutually exclusive exclusion counts.
labels=['Case Reports','Editorial','Letter','Comment']
write('04_出版类型多标签与互斥排除对照.csv',[{'publication_type':t,'all_records_with_tag':sum(t in r['publication_types'] for r in rs),'legacy_excluded_with_tag':sum(t in r['publication_types'] and r['decision']=='exclude' for r in rs),'legacy_included_with_tag':sum(t in r['publication_types'] and r['decision']=='include' for r in rs),'not_additive':'标签可能重叠，不能相加作为总排除数'} for t in labels],['publication_type','all_records_with_tag','legacy_excluded_with_tag','legacy_included_with_tag','not_additive'])
write('05_年度缺项筛选影响.csv',[{'year':y,'legacy_eligible':sum(r['year']==y for r in eligible),'field_complete':sum(r['year']==y for r in complete),'note':'期刊年份字段；并非最终首次发表年份统计'} for y in sorted(set(r['year'] for r in eligible))],['year','legacy_eligible','field_complete','note'])
summary={'status':'PRE_SUBMISSION_AUDIT_NOT_FINAL_RESULTS','audited_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'raw_pubmed':len(rs),'original_query':sum(r['retrieval_group']=='original' for r in rs),'expanded_additions':sum(r['retrieval_group']=='additional' for r in rs),'legacy_decisions':dict(Counter(r['decision'] for r in rs)),'legacy_eligible':len(eligible),'missing_counts_in_legacy_eligible':dict(Counter(k for r in eligible for k in missing(r))),'legacy_eligible_with_any_core_missing':sum(bool(missing(r)) for r in eligible),'field_complete_legacy_eligible':len(complete),'field_complete_original_query_legacy_eligible':sum(r['retrieval_group']=='original' for r in complete),'field_complete_expansion_additions':sum(r['retrieval_group']=='additional' for r in complete),'legacy_eligible_with_reference_lists':sum(bool(r.get('references')) for r in eligible),'wos_files':len(list((R/'input/wos').glob('*.txt'))),'mutually_exclusive_legacy_exclusions':dict(reasons)}
assert len(rs)==sum(summary['legacy_decisions'].values())
assert sum(reasons.values())==summary['legacy_decisions']['exclude']
assert len(eligible)==len(complete)+summary['legacy_eligible_with_any_core_missing']
assert len(rs)==len({r['pmid'] for r in rs})
(O/'audit_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
(O/'06_本轮核验结论.md').write_text(f'''# SCI投稿前数据核验结论

本轮核验的是已保存的604条PubMed记录，没有重新执行数据库检索，也没有取得WoS数据；因此不将本文件称为双库筛选结果。

旧筛选纳入广义范围文献173篇，其中{summary['legacy_eligible_with_any_core_missing']}篇至少缺少一项本次要求的核心书目字段；关键词缺失51篇、机构缺失4篇，存在重叠。现有字段完整者为{len(complete)}篇，其中原式来源{summary['field_complete_original_query_legacy_eligible']}篇、扩展新增来源{summary['field_complete_expansion_additions']}篇。这些数量仅为补全前字段核验，不是本研究最终纳入数。应先跨库补全再判定排除，不能直接把缺项记录全部删掉并宣布终稿完成。

旧纳入173篇中仅{summary['legacy_eligible_with_reference_lists']}篇保存有参考文献列表，而且“有列表”不表示参考文献完整。现阶段不能将其用于全样本共被引排名或期刊连接强度推断。

旧604条的PMID均唯一，且纳入、明确排除、信息不足三者之和等于604；415条互斥排除原因计数核对通过。病例报告、Editorial及Letter多标签计数另表报告，避免把类型标签与最终排除原因混加。

## 本轮已完成

逐条输出604条核心字段核验；列出旧纳入但缺项的文献；区分互斥排除与类型多标签；生成年度缺项影响表；记录源文件SHA-256及复现脚本。

## 尚不能宣布“可投稿终稿”的原因

WoS文件为0；未能建立数据库浏览器连接，公开页面返回403。CiteSpace尚无可核验的正式运行工程。缺少双库去重匹配、跨库补全、新标准最终筛选及正式网络结果。正文可以润色，但以上数据条件不能由语言修订替代。

接续只需优先提供WoS Core Collection全部批次Plain Text（Full Record and Cited References），同时给出实际执行检索式、命中数、索引与检索日期。作者、单位、基金和利益冲突等署名信息可在分析完成后统一填入，不妨碍当前数据准备。
''',encoding='utf-8')
print(json.dumps(summary,ensure_ascii=False,indent=2))
