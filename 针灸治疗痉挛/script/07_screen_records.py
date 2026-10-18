"""Apply the recorded title/abstract decisions, not a keyword-based inclusion classifier.

Index lists below refer to frozen PMID review-order files. Every non-rule-based
decision is exported by PMID with its supporting text and rationale.
"""
from pathlib import Path
import json,re,csv,collections
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'output/06_完整分析'
records=json.loads((OUT/'all_records.json').read_text());byid={r['pmid']:r for r in records}
decisions={}
original_exclusions={
0:('Study scope','脑瘫中医干预的广泛文献计量研究，未提供针灸抗痉挛的独立临床或机制分析。'),
3:('Intervention','比较冲击波与穴位经皮电刺激，无穿刺针灸或艾灸干预。'),
7:('Topic','主要研究神经根性疼痛及炎症，痉挛仅为疼痛相关步态描述，未明确评价抗痉挛效应。'),
17:('Intervention','组间比较rTMS；针灸仅为两组共同的背景康复措施，未评价针灸。'),
26:('Intervention','仅研究经皮穴位电刺激。'),
35:('Intervention','研究经皮神经电刺激，非穿刺针灸或艾灸。'),
36:('Intervention','仅研究经皮穴位电刺激。'),
38:('Study scope','非药物疗法整体文献计量研究，未独立评价针灸抗痉挛临床或机制问题。'),
42:('Topic','整体卒中神经功能综述；摘要只列举痉挛，未报告针对痉挛的独立针灸分析。'),
47:('Topic','开放原文核查：针灸出现在检索词和讨论背景，纳入干预未形成针灸抗痉挛分析。'),
55:('Topic','经筋筋膜理论展望，痉挛为理论性选穴描述。'),
56:('Intervention','振动加热砭石与运动干预，不属于本次针刺、干针或艾灸范围。'),
61:('Intervention','耳穴王不留行籽按压，无针刺或艾灸。'),
63:('Topic','一般卒中神经可塑性综述；痉挛仅为背景中的后遗症，未明确独立抗痉挛分析。'),
78:('Topic','评价卒中步态，未明确痉挛人群或痉挛结局。'),
81:('Intervention','激光穴位刺激，无穿刺针灸或艾灸。'),
82:('Protocol','摘要使用we will及计划方法，未报告完成的试验结果。'),
85:('Protocol','摘要明确计划招募316例，结局以将来时描述，无完成结果。'),
87:('Topic','一般电针不良事件综述，痉挛为不良事件而非治疗对象。'),
90:('Topic','评价减压手术疗效；针灸为术前既往治疗，非所评估治疗。'),
92:('Intervention','针刀微创松解，不属于本次干预范围。'),
95:('Topic','牙拔除术后张口受限，未明确痉挛诊断或肌痉挛结局。'),
97:('Topic','评价微血管减压预后，既往针灸仅为预测变量。'),
99:('Topic','Oddi括约肌功能障碍综述，针灸针对疼痛，痉挛为药物不良反应背景。'),
101:('Case report','摘要为单名患者的干针治疗，数据库未标出Case Reports。'),
106:('Topic','一般缺血性卒中基础机制综述，抗痉挛仅在引言背景提及。'),
109:('Topic','研究多发性硬化步行，痉挛只是既往文献背景，未测量痉挛。'),
113:('Intervention','激光穴位刺激。'),
114:('Topic','急性腰扭伤的肌肉力位移与温度研究，无明确痉挛诊断或抗痉挛分析。'),
117:('Intervention','仅研究经皮穴位电刺激。'),
119:('Topic','偏瘫fMRI研究，无明确痉挛人群或痉挛测量，痉挛为结论外推。'),
122:('Intervention','穴位埋线，不属于本次针刺、干针或艾灸范围。'),
129:('Topic','肌筋膜疼痛干预综述；痉挛在激光治疗假说中出现。'),
132:('Intervention','综述纳入标准明确排除针灸及干针。'),
133:('Topic','贝尔麻痹物理治疗，未提供针灸针对面肌痉挛的可区分结果。'),
134:('Topic','下颌脱位复位研究，肌痉挛仅为辅助镇静镇痛的解释，无抗痉挛结果。'),
135:('Topic','前列腺炎/盆腔疼痛综合征分型综述，未单独分析针灸治疗盆底痉挛。'),
136:('Intervention','仅研究经皮穴位电刺激。'),
140:('Population','健康受试者脑响应实验，面肌痉挛仅为疗效推测。'),
142:('Population','健康志愿者肌电研究，无痉挛患者或病理模型。'),
144:('Topic','肌筋膜疼痛试验，痉挛仅为疾病背景描述。'),
146:('Intervention','电极穴位刺激与按压，无针刺或艾灸。'),
150:('Topic','肌筋膜触发点综述，痉挛治疗指触发点注射，非所评估干预。'),
151:('Study scope','治疗手册制定方法文章，未报告已完成痉挛治疗或机制研究。'),
157:('Intervention','评价鞘内κ受体激动剂，电刺激作用是与既往研究的机制比较。'),
158:('Intervention','皮肤电极HANS刺激，无穿刺针灸或艾灸。'),
160:('Intervention','放疗后电疗与高压氧，非针灸抗痉挛研究。'),
161:('Topic','恶性肿瘤疼痛主题，题录未显示直接抗痉挛研究。'),
162:('Unavailable text','题名相关，但PubMed无摘要；公开检索未取得可核实摘要或全文。另列题名相关记录敏感性清单。'),
163:('Population','志愿者反射实验，未建立痉挛患者或疾病模型。'),
164:('Topic','血液流变学理论讨论，痉挛只作针灸可能作用的示例。'),
165:('Topic','弹响指康复，非痉挛疾病研究。'),
166:('Topic','节段性疼痛综述，未明确独立抗痉挛研究。'),
168:('Topic','混合疼痛门诊病例系列，摘要未识别痉挛亚组或抗痉挛结局。')}
additional_exclusions={
9:('Intervention','超声引导针刀，非本次针刺、干针或艾灸范围。'),
16:('Case report','题名明确单病例报告。'),
17:('Topic','颈背肌筋膜疼痛试验，痉挛为背景疾病特征，未报告痉挛评价。'),
19:('Population','无神经或骨科损伤的潜在触发点受试者，无痉挛疾病。'),
21:('Topic','痛经动物操作方法；spasmodic描述腹痛，无单独痉挛病理或功能评价。'),
22:('Protocol','摘要明确will be conducted/随机分组计划，无完成结果。'),
25:('Topic','急性腰痛指南；痉挛对应肌松药，非针灸抗痉挛评价。'),
26:('Topic','神经根性疼痛机制综述，痉挛仅为镇痛假说，未独立评价。'),
38:('Population','离体正常小鼠肌肉安全性实验，无痉挛疾病模型。'),
39:('Population','正常尸体针刺安全解剖研究，无痉挛病理模型或患者。'),
41:('Intervention','经络按摩综述，非针刺、干针或艾灸。'),
43:('Topic','肠易激综合征治疗综述；spasmolytic为药物类别，非针灸抗痉挛评价。'),
44:('Population','正常尸体解剖安全研究，无痉挛患者或疾病模型。'),
46:('Protocol','儿童按摩研究方案，明确未来招募和计划结局。'),
49:('Topic','HTLV血液供者流行病学，针灸仅为感染危险因素。'),
53:('Protocol','摘要明确未来随机分组及计划结局，未报告结果。'),
61:('Population','健康志愿者食管测压，痉挛性疾病为未来应用推测。'),
63:('Topic','便秘与IBS指南，未独立评价针灸抗痉挛研究。'),
64:('Study scope','近期肌筋膜文献评述，未提供针灸/干针对痉挛的独立分析。'),
66:('Topic','肉毒毒素腰痛研究，针灸为比较措施，未评价痉挛病理或结局。'),
67:('Intervention','硫代秋水仙碱研究，干针仅为背景疗法。'),
68:('Topic','神经病理性疼痛治疗，spasmodic描述疼痛性质，不是肌痉挛。'),
69:('Unavailable text','题名相关；PubMed无摘要，公开检索未取得可核实摘要或全文，列敏感性清单。'),
72:('Topic','多发性硬化生活质量试验，肌痉挛作为不良反应而非治疗对象。'),
74:('Unavailable text','30例治疗题名相关，非仅凭病例系列排除；未取得可核实摘要或全文，列敏感性清单。'),
75:('Intervention','指压，无针刺或艾灸。'),
76:('Case report','题名及摘要明确单病例报告。'),
77:('Topic','胃肠病补充疗法综述；抗痉挛描述对应薄荷油，针灸对应恶心呕吐。'),
79:('Protocol','肉毒毒素试验方案，未报告完成研究结果。'),
81:('Intervention','按压与骨科手法，无针刺或艾灸。'),
82:('Topic','高血压模型与微循环研究，血管痉挛仅为形态背景，不为独立治疗问题。'),
83:('Topic','混合疼痛诊疗报告，无法区分针灸治疗痉挛性斜颈的结果。'),
85:('Topic','穴位针刺诱发痉挛的诊断现象，非抗痉挛治疗。'),
87:('Unavailable text','无PubMed摘要；公开记录为广泛中枢疾病综述，未取得可核查的目标疾病全文。')}
for group,exclusions in [('original',original_exclusions),('additional',additional_exclusions)]:
    ids=json.loads((OUT/f'{group}_review_order.json').read_text())
    for i,pmid in enumerate(ids):
        if i in exclusions:
            reason,note=exclusions[i];decisions[pmid]={'decision':'exclude','reason':reason,'note':note}
        else:
            decisions[pmid]={'decision':'include','reason':'Eligible','note':'题名/摘要核查：直接研究针刺、干针或艾灸治疗痉挛/痉挛状态，或对该问题开展证据综述、机制/治疗使用研究。'}
prior=json.loads((ROOT/'script/relevance_decisions.json').read_text())
for pmid,d in prior.items():
    if d['decision']=='exclude':decisions[pmid]={'decision':'exclude','reason':'Topic','note':d['reason']}
manual_types={'4414100':('Letter','题名核实信函。'),'38188810':('Retraction notice','题名核实撤稿通知。'),'30377465':('Case report','摘要核实两例病例报告。')}
for pmid,(reason,note) in manual_types.items():decisions[pmid]={'decision':'exclude','reason':reason,'note':note}
for r in records:
    pmid=r['pmid']; t=r['title']; types=r['publication_types']
    if r['primary_exclusion']:
        d={'decision':'exclude','reason':r['primary_exclusion'],'note':'依据PubMed出版类型或发表语言字段。'}
    elif any(x['type']=='RetractionIn' for x in r.get('comments_corrections',[])):
        d={'decision':'exclude','reason':'Retracted Publication','note':'PubMed关联撤稿通知。'}
    elif re.search(r'\bprotocol\b',t,re.I):
        d={'decision':'exclude','reason':'Protocol','note':'题名明确为研究方案，无已完成研究结果。'}
    elif pmid in decisions:d=decisions[pmid]
    else:raise ValueError('Unreviewed: '+pmid+' '+t)
    if d['reason'] in ['Case Reports','Case report']:d['reason']='Case report'
    if d['reason']=='Book chapter':d['reason']='Book/report'
    r['decision']=d['decision'];r['exclusion_reason']='' if d['decision']=='include' else d['reason'];r['decision_note']=d['note']
    r['review_method']='AI-assisted title/abstract review; no independent human dual screening'
    sentences=re.split(r'(?<=[.!?])\s+|\n',r['abstract'])
    r['supporting_text']=' '.join(s for s in sentences if re.search('spas|ashworth|tardieu|muscle tone|will be',s,re.I)) or r['title']
    r['supporting_source']=r['url']
    if pmid=='38510783':r['supporting_source']='https://europepmc.org/articles/PMC10953734'
(OUT/'screened_records.json').write_text(json.dumps(records,ensure_ascii=False,indent=2))
fields=['pmid','doi','title','year','journal','languages','publication_types','retrieval_group','decision','exclusion_reason','decision_note','supporting_source','supporting_text','review_method']
def output(name,rows):
    with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fields);w.writeheader();w.writerows({k:json.dumps(r[k],ensure_ascii=False) if isinstance(r[k],list) else r[k] for k in fields} for r in rows)
output('01_604条逐篇筛选.csv',records);output('02_纳入文献.csv',[r for r in records if r['decision']=='include'])
output('03_未取得摘要全文记录.csv',[r for r in records if r['exclusion_reason']=='Unavailable text'])
print(collections.Counter(r['decision'] for r in records));print(collections.Counter(r['exclusion_reason'] for r in records if r['decision']=='exclude'))
