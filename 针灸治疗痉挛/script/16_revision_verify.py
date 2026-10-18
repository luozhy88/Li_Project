"""Verify the revised analysis and write reproducible audit outputs."""
from pathlib import Path
import json,csv,hashlib,collections,re,zipfile
ROOT=Path(__file__).resolve().parents[1];B=ROOT/'output/09_审稿后重分析';O=ROOT/'output/10_审稿后修订稿'
r=json.loads((B/'screened_records.json').read_text());main=json.loads((B/'main/analysis_records.json').read_text());broad=json.loads((B/'broad/analysis_records.json').read_text());s=json.loads((B/'main/analysis_summary.json').read_text());b=json.loads((B/'broad/analysis_summary.json').read_text());by={x['pmid']:x for x in r}
assert len(r)==len({x['pmid'] for x in r})==604
assert len(main)==146 and len(broad)==173
assert sum(x['decision']=='exclude' for x in r)==415
assert sum(x['decision']=='unassessable' for x in r)==16
assert 415+16+len(broad)==604 and len(broad)-len(main)==27
assert all(by[k]['decision']=='include' for k in ['34777214','37667737','34219954','25780676'])
assert all(by[k]['reviewed_design']=='Animal experiment' for k in ['42022268','41503987','41299887'])
assert by['32815384']['reviewed_design']=='Other human research'
assert by['42037903']['reviewed_design']=='Other human research'
assert by['38013386']['reviewed_design']=='Uncertain clinical design'
assert sum(s['designs'].values())==sum(s['diseases'].values())==sum(s['annual'].values())==146
assert s['keyword_papers']==sum(bool(x['terms']) for x in main)==109
assert s['original_included']+s['additional_included']==146
assert s['excluded']==415 and s['information_insufficient']==16 and s['eligible_outside_analysis']==27
assert s['designs']==dict(collections.Counter(x['design'] for x in main))
for sub,summary in [('main',s),('broad',b)]:
 d=B/sub;cross=list(csv.DictReader((d/'12_干预专题设计疾病交叉表.csv').open(encoding='utf-8-sig')))
 assert sum(int(x['papers']) for x in cross)==summary['included']
 for x in cross:assert sum(int(x[k]) for k in summary['designs'])==int(x['papers'])
 assert len(list(csv.DictReader((d/'network_terms_nodes.csv').open(encoding='utf-8-sig'))))==summary['network_stats']['terms']['nodes']
 assert len(list(csv.DictReader((d/'network_terms_edges.csv').open(encoding='utf-8-sig'))))==summary['network_stats']['terms']['edges']
 for file in d.glob('figures/*.png'):assert file.stat().st_size>5000
components=json.loads((O/'manuscript_components.json').read_text());body=components['body']
assert [int(i) for i in re.findall(r'\[FIG(\d+)\]',body)]==list(range(1,12))
assert [int(i) for i in re.findall(r'\[TABLE(\d+)\]',body)]==list(range(1,8))
assert len(components['figures'])==11 and len(components['tables'])==7
report={'status':'passed','unique_records':604,'excluded':415,'information_insufficient':16,'broad':173,'main':146,'scope_only':27,'keywords':109,'design_counts':s['designs'],'checks':['record accounting','corrected type conflicts','study classification totals','cross-tabulation margins','network node/edge counts','all figure and table anchors sequential'],'limitations':'Numerical consistency is not independent human eligibility review.'}
(B/'validation_report.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
query=json.loads((ROOT/'output/01_检索与研究方案/raw/pubmed_sensitivity_esearch.json').read_text())['esearchresult']['querytranslation']
supp=f'''# 补充方法与复现说明

修订日期：2026-10-18。主文所有图表位于正文。

## 数据队列与分母

固定检索记录604条；明确排除415条；信息不足16条；广义符合标准173篇，其中神经性痉挛主分析146篇，其他痉挛领域27篇。信息不足不等于临床不合格，未对其研究设计或结局进行猜测。

主分析：原式97篇，新增49篇；关键词109篇；机构字段144篇；国家/地区可定位143篇；机构体系139篇。广义集合：关键词122篇。全计数允许主体贡献之和超过文献总数。所有数量为出版物数，不是独立试验数。

## 完整数据库检索转换式

```text
{query}
```

原式与扩展式实际返回结果保存在output/01_检索与研究方案/raw/；604个固定PMID及原始XML在原有数据包中保留。检索保存日期2026-10-11，发表截止2026-10-10，新增209条XML补取日期2026-10-18。没有把后续实时新增命中混入队列。

## 规则修改及证据

主分析限定神经性痉挛，其他疾病保留为范围敏感性集合；这是修订阶段明确的范围，不声称预注册。保留针刺作为明确研究治疗方案的联合治疗和前后观察，排除只在常规背景治疗清单中出现的针灸。干针专题按题名dry needl词形识别，并纳入开放全文证实的PMID 35024159针介导肌内电刺激；其余专题可能包含干针讨论，不能当成互斥治疗组。

补入PMID 34777214、37667737、34219954及25780676，最后一篇仅进入广义集合。修正PMID 42022268、41503987、41299887为动物实验；42037903和32815384为其他人体研究；38013386保留设计冲突。具体来源摘要见03_604条筛选与结构化依据.csv，原判断与新判断见01_审稿修改日志.csv。病例报告排除标签复核清单共91条。支持文本与原始摘要用于证据定位，不替代独立审稿人。

信息不足清单见04_信息不足记录.csv。其中部分为无摘要旧文，另有摘要不足以核实统一设计的文章。并未声称穷尽付费全文渠道或联系作者。获得可核实原文后，应更新筛选依据并重新计算，不只修改主文数量。

## 网络与分类敏感性

主分析、原式子集、非综述子集、卒中子集均在频次2、3、5下计算；广义集合另行计算。删除设计词的名单为meta-analysis、systematic review、network meta-analysis、randomized controlled trial、review、clinical trial、pilot study、case series、randomized trial、scoping review。关键词分母保持109篇。

PMID 38013386归入随机临床研究时，该类43→44篇；归入其他人体研究时，该类19→20篇。其他计数与关键词网络不变。未使用自定义突现模型，也未进行年度线性拟合、未来预测或时期差异显著性检验。

国家补全采用明确机构白名单及地址地名；逐条映射见main/05_机构与国家映射.csv。关键词规范化见main/04_关键词规范化词典.csv。作者署名未全面身份消歧。软件为Python 3.11.7、NetworkX 3.1、NumPy 1.26.4、SciPy 1.11.4、Matplotlib 3.8.0；图为PNG与SVG，随机种子42。

## 复现顺序

```sh
python3 script/12_revision_screen.py
/opt/anaconda3/bin/python script/13_revision_analysis.py --broad
/opt/anaconda3/bin/python script/13_revision_analysis.py
python3 script/14_revision_manuscript.py
/Users/apple/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 script/15_revision_docx.py
python3 script/16_revision_verify.py
```

12脚本读取原版screened_records.json与analysis_records.json并应用可追溯修订，旧文件不被覆盖。13重新生成主分析与广义集合。14从重算结果生成文章组件，15将7张表及11幅图依次嵌入正文。文档排版工具依赖python-docx及已安装的documents技能辅助库；其他机器需要配置等价依赖和路径。

## 研究者责任

本轮为AI辅助复核，未完成独立人类双筛或完整风险偏倚评价。不得将自动验证通过写作临床资格已获得独立验证。作者、单位、基金、利益冲突、贡献和伦理适用性需由署名研究者确认；当前文件不虚构这些声明。
'''
(O/'03_补充方法与复现说明.md').write_text(supp)
print(json.dumps(report,ensure_ascii=False))
