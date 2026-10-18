"""Versioned reviewer corrections; source dataset is preserved."""
from pathlib import Path
import json,csv,re,collections
ROOT=Path(__file__).resolve().parents[1];OLD=ROOT/'output/06_完整分析';OUT=ROOT/'output/09_审稿后重分析';OUT.mkdir(exist_ok=True)
rs=json.loads((OLD/'screened_records.json').read_text());prior={r['pmid']:r for r in json.loads((OLD/'analysis_records.json').read_text())}
include={
 '34777214':('Randomized clinical study','Stroke','210名卒中患者分为干针、假干针和对照三组；摘要报告随机分配，按实际内容覆盖病例报告标签。'),
 '37667737':('Other human research','Cerebral palsy','痉挛型脑瘫83名接受针灸推拿、85名接受常规康复；摘要未说明随机，归为非随机/未报告随机的人体研究。'),
 '34219954':('Other human research','Multiple sclerosis','12名多发性硬化患者接受12次干针，治疗前后评价PSFS及功能；系统性前后观察，非单病例。'),
 '25780676':('Other human research','Facial spasm/synkinesis','32名面肌痉挛患者接受针刺联合胎盘药物穴位注射，有统一量表及前后统计；仅纳入广义范围，不能分离针刺效应。')}
fixdesign={'32815384':'Other human research','42022268':'Animal experiment','41503987':'Animal experiment','41299887':'Animal experiment','42037903':'Other human research','38013386':'Uncertain clinical design'}
# A neurological-spasticity population, not merely any paper containing spasm.
neuro_diseases={'Stroke','Cerebral palsy','Multiple sclerosis','Brain injury/disorders of consciousness','Spinal cord disorders','Mixed neurological disorders'}
neuro_mixed={'35264894','30408516','29042299','42062015','33399051'}
# No available abstract: do not infer an observational design from the Case Reports tag.
unavailable={'18246691','12400435','11789328','8246597','1598002','1861516','3926406','6567732','7274598','15510797','2866668','6145300'}
intervention_exclude={'40098901':'穴位冲击波，不属于穿刺针刺或艾灸。','6416756':'题名为电疗与法拉第刺激，未证实穿刺针刺。','2576564':'穴位按摩，不属于穿刺针刺或艾灸。','30037913':'金线植入，属于植入材料干预而非本研究普通针刺。','1978504':'研究气功能量存储，针灸仅为背景比较。','2568074':'研究气功状态，非针刺干预。'}
changes=[];case_audit=[]
for r in rs:
 old=r['decision'];oldreason=r['exclusion_reason'];id=r['pmid'];r['old_decision']=old;r['old_reason']=oldreason
 if id in include:
  design,disease,note=include[id];r.update(decision='include',exclusion_reason='',decision_note=note)
 elif id in intervention_exclude:r.update(decision='exclude',exclusion_reason='Intervention',decision_note=intervention_exclude[id])
 elif id in unavailable:r.update(decision='unassessable',exclusion_reason='Unavailable text',decision_note='现有题录或摘要不足以核实设计及目标干预；不因病例报告标签直接判定不合格。')
 elif oldreason=='Unavailable text':r['decision']='unassessable'
 if 'eng' not in r['languages']:
  r.update(decision='exclude',exclusion_reason='Non-English',decision_note='按全文发表语言字段排除，英文摘要不等于英文论文。')
 if r['decision']=='include':
  if id not in include:design,disease=prior[id]['design'],prior[id]['disease']
  design=fixdesign.get(id,design)
  r['reviewed_design']=design;r['reviewed_disease']=disease
  r['main_scope']=disease in neuro_diseases or id in neuro_mixed
  # Explicitly broad mechanism reviews are scope-sensitivity records, not primary neurological papers.
  r['target_domain']='Neurological spasticity' if r['main_scope'] else ('Visceral/vascular spasm' if disease=='Visceral/vascular spasm' else 'Other/mixed spasm disorders')
  r['underlying_disease']='Stroke' if id=='34874522' else disease
  r['dry_needling_focus']=bool(re.search('dry needl',r['title'],re.I)) or id=='35024159'
  r['focus_basis']='title' if id!='35024159' else 'open full text: needle-based intramuscular stimulation'
  r['classification_evidence']=r['abstract']
  r['decision_note']=include[id][2] if id in include else f"研究对象：{disease}；设计：{design}；范围：{r['target_domain']}。资格依据为本条摘要及所列支持文本；不据此评价疗效质量。"
  if id=='29436374':r['decision_note']+=' 针灸单用与针灸联合贴扎均为主动评价的方案，纳入出版物地图；不解释为针灸独立疗效证据。'
  if id=='38013386':r['decision_note']+=' 题名称观察性，摘要/索引称随机；保留设计不确定，不任意归类。'
 else:r['main_scope']=False
 r['review_method']='AI-assisted reviewer revision; no independent human dual screening'
 if oldreason=='Case report':case_audit.append({'pmid':id,'title':r['title'],'old_reason':oldreason,'new_status':r['decision'],'new_reason':r['exclusion_reason'],'evidence':r['abstract'] or r['title'],'note':r['decision_note']})
 if old!=r['decision'] or oldreason!=r['exclusion_reason'] or id in fixdesign:changes.append({'pmid':id,'title':r['title'],'old_decision':old,'old_reason':oldreason,'new_decision':r['decision'],'new_reason':r['exclusion_reason'],'note':r['decision_note']})
(OUT/'screened_records.json').write_text(json.dumps(rs,ensure_ascii=False,indent=2))
def save(name,rows):
 with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:
  w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)
save('01_审稿修改日志.csv',changes);save('02_91条病例标签复核.csv',case_audit)
fields=['pmid','title','decision','exclusion_reason','main_scope','reviewed_design','reviewed_disease','target_domain','underlying_disease','dry_needling_focus','focus_basis','decision_note','supporting_text','abstract','url','review_method']
save('03_604条筛选与结构化依据.csv',[{k:r.get(k,'') for k in fields} for r in rs])
save('04_信息不足记录.csv',[{k:r.get(k,'') for k in fields} for r in rs if r['decision']=='unassessable'])
print(collections.Counter(r['decision'] for r in rs));print('main',sum(r['main_scope'] for r in rs));print('case reviewed',len(case_audit))
