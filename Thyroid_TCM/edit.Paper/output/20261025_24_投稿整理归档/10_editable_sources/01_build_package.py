from pathlib import Path
import re, shutil, subprocess, json
from docx import Document
from docx.shared import Inches,Pt
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
root=Path('/home/data1/zhiyu/data/Github/Li_Project/Thyroid_TCM/edit.Paper/output')
out=root/'20261025_19_JCPSP_submission'; src=out/'10_editable_sources'
s=(root/'20261025_21_article_version_reconciled.md').read_text(); title=s.splitlines()[0][2:]
auth=s.split('# Abstract')[0]
main='# '+title+'\n\n# Abstract'+s.split('# Abstract',1)[1]
# Keep identity details separate from the review manuscript.
main=main.replace('This study was conducted at Shanghai University of Medicine & Health Sciences, Shanghai, China. ','')
main=main.replace('**Methods:**','**Study design:** Secondary analysis of public databases and paired tumour–normal expression data.\n\n**Methodology:**',1)
# Place/duration is pending author information, documented in the submission guide, not invented.
captions={int(n):t for n,t in re.findall(r'\*\*Figure (\d+):\*\* ([^\n]+)',main)}
main=re.sub(r'!\[\]\([^\n]+\n','',main)
main=re.sub(r'\*\*Figure \d+:\*\* [^\n]+\n','',main)
mapping={1:'Supplementary Figure S1',2:'Figure 1',3:'Figure 2',4:'Figure 3',5:'Supplementary Figure S2'}
main=re.sub(r'Figure ([1-5])',lambda m:mapping[int(m[1])],main)
# Declarations on title page; non-identifying ethics/data statements retained.
declarations=[]
for heading in ['Conflict of Interest','Funding','Author Contributions']:
 pat=r'# '+heading+r'\n\n(.*?)(?=\n# |\Z)'
 match=re.search(pat,main,re.S)
 if match: declarations.append(match[0]); main=main[:match.start()]+main[match.end():]
main+='\n# Figure legends\n\n'+'\n\n'.join('**Figure '+str(i)+':** '+re.sub(r'Figure ([1-5])',lambda m:mapping[int(m[1])],captions[old]) for i,old in [(1,2),(2,3),(3,4)])+'\n'
titlepage=auth+'\nArticle type: Original Article (proposed; suitability subject to editorial assessment).\n\nRunning title: Xiao Ying Tang core herbs and papillary thyroid cancer\n\n'+ '\n\n'.join(declarations)
cover=f'''# Cover letter — author review draft

Dear Editor,

Please consider our manuscript, “{title}”, for publication as an Original Article in the Journal of the College of Physicians and Surgeons Pakistan.

The study examines genes associated with core herbs selected from Xiao Ying Tang and evaluates their disease relevance in an independent papillary thyroid carcinoma expression dataset. The analysis uses matched tumour and normal thyroid tissue from 24 patients, with repeated tissue blocks averaged within each patient and tissue type. Candidate-set enrichment is assessed against both the full measured-gene background and a disease-restricted background.

XYT- and Xia Ku Cao-associated genes were over-represented among differentially expressed genes under both backgrounds. The smaller Zhe Bei Mu set was not significantly enriched, and none of the three preselected genes met the combined differential-expression criteria. We report these findings as disease associations that guide functional investigation, without claiming treatment efficacy or herb synergy.

The manuscript includes three main figures, with supplementary gene lists, expression results, sample manifests, enrichment results and supporting illustrations. The work uses public, de-identified data and involved no new participant recruitment or intervention. All authors have reviewed and approved the manuscript. Funding, author contributions and conflicts of interest are disclosed on the title page.

Thank you for considering this work.

Yours sincerely,

Xing Zhang\n\nCorresponding author\n\nDepartment of Respiratory Disease, Shuguang Hospital Affiliated to Shanghai University of Traditional Chinese Medicine, Shanghai 201203, China\n\nEmail: zhangxing2012@shutcm.edu.cn\n\nTelephone: 021-20256375
'''
inquiry=f'''# Pre-submission enquiry — not sent

To: publications@cpsp.edu.pk

Subject: Pre-submission enquiry: computational thyroid cancer study and equal-contribution authorship

Dear Editorial Office,

We are considering submitting “{title}” to JCPSP. The study combines network pharmacology with a secondary analysis of an independent public transcriptome dataset comprising paired papillary thyroid carcinoma and normal thyroid tissue from 24 patients. It includes no new participant recruitment, intervention, animal experiments or wet-laboratory experiments.

Before formal submission, could you please advise on the following points?

1. Would this design be eligible for consideration as an Original Article?
2. Our five-author manuscript identifies Tao Li, Jiacai Li and Yajie Ding as equal first contributors, with Xing Zhang as the sole corresponding author. Would an equal-contribution statement be permitted under your current authorship policy?
3. For this secondary analysis of public, de-identified data, what institutional exemption documentation would your office require?
4. The underlying public dataset was published in 2015; the present analysis is a separate study. How should the study-duration and data-age requirements be applied in this case?

5. AI assistance has been used during manuscript editing and review. We have noted your policy prohibiting AI-assisted editing and polishing. Does this make the manuscript ineligible for consideration? We do not wish to submit it in conflict with your policy.

We would be grateful for your guidance before preparing the signed submission forms.

Yours sincerely,

Xing Zhang\n\nCorresponding author\n\nEmail: zhangxing2012@shutcm.edu.cn
'''
guide='''# JCPSP 投稿步骤与核对清单

整理日期：2026年10月25日。当前状态：投稿准备包，尚未投稿，不建议直接投 JCPSP。当前官网禁止使用 AI 生成、编辑和润色文本，与本文实际修改过程冲突；另有共同第一作者限制。不能将本文件夹视为可直接提交的合规终稿。

## 一、登录哪个网站

期刊官网：https://www.jcpsp.pk/

官方作者指南：https://www.jcpsp.pk/instructions-to-author.php

官方指南指定投稿入口：https://www.jcpsp.pk/oas

编辑部联系邮箱：publications@cpsp.edu.pk

从期刊官网进入投稿系统。不要误选 JCPSP Postgraduate 或 JCPSP Case Reports；它们不是本次目标期刊。核查时投稿系统页面未能加载，因此下列为操作顺序说明，并非已登录账户实测的按钮教程；界面字段以登录后显示为准。

## 二、先解决的事项

官网当前禁止 AI 或数字语言工具生成、编辑和润色作者文本。本稿已接受 AI 辅助修改，删除声明不改变这一事实。建议改选允许如实披露辅助工具的期刊；若坚持 JCPSP，应先向编辑部如实说明并取得政策解释。

官网当前限定单一第一作者，与本稿三位共同第一作者冲突。文件保留现有署名，没有擅自删除共一。先用 08 英文预询问信征询编辑部；本信未发送。

官网要求研究文章提交机构伦理批准或豁免文件，附编号及日期。本文没有新受试者研究，但不能替机构签发豁免。请由研究单位提供适用文件，或取得编辑部对替代说明的书面认可。

作者指南对作者的医学/博士及相关专业资质有要求。五位作者的学位和专业资格需逐一确认；现有资料不能替代资格审核。仅有数据整理或指导的贡献描述，也需由作者核对是否覆盖期刊要求的完整作者责任。

摘要要求研究设计、地点和起止时间。研究单位已明确；真实研究起止年月尚未提供。本版没有杜撰日期，正式上传前需要补齐相应摘要栏目。GSE60542 数据来源年份不能充当本研究开展日期。

现有24条文献保留了方法原始出处。近五年文献占比尚未达到指南要求的多数，正式投稿前需按最终文章论点更新相关临床及机制文献，保留实际所用软件和数据集的原始引用。不能用新版数据库文献暗示实际使用了新版数据。

作者尚需确认未发表、未一稿多投及任何预印本状态，并在官方 ACP 表上签名。现有介绍信没有代作者作出尚未确认的原创性/独家投稿承诺。全体作者批准最终稿已确认。

## 三、文件怎么用

01 英文正文：用于同行评审，作者身份页单列；主文保留3幅图的引用及图注，图片见09文件夹。正文保留所有关键统计结果，阴性结果不删除。

02 作者信息页：姓名、顺序、共同第一作者、单位、通讯作者、资助、利益冲突及贡献。正式提交前补齐其他作者的联系信息。

03 英文介绍信：给编辑的 cover letter，需作者审阅后使用。

04 中文审阅稿：最新完整中文对照稿，仅供团队审阅；仍沿用原5图编号，不作为英文投稿附件上传。英文稿已调整为3幅主图及2幅补充图，内容相同，编号对应见补充图文档。

05 补充材料：结果表、基因列表、样本信息及质量控制等，不含分析代码。06 补充图：S1基因交集图、S2探索性对接图及图注。

07 本文档：团队操作说明，不上传为论文附件。

08 预询问信：用于解决范围、署名及伦理文件要求；不等同于正式投稿。

09 figures：主图1 GO、主图2 KEGG、主图3独立转录组。原始图像按现有分辨率保留；如系统要求 TIFF/EPS 或特定 dpi，提交前按源图实际尺寸导出，不能仅改扩展名或伪造分辨率。

10 editable sources：Markdown 可编辑来源，仅供内部修改，不上传。

## 四、正式投稿步骤

1. 通讯作者 Xing Zhang 使用本人邮箱 zhangxing2012@shutcm.edu.cn，从官网进入官方投稿系统。若无账户，按页面提示注册并验证邮箱；已有账户直接登录。官网要求由通讯作者本人提交。
2. 开始新投稿；在编辑部确认本文类型后选择 Original Article。填写英文标题、摘要和关键词，内容与最终 Word 一致。
3. 按 Tao Li、Jiacai Li、Yajie Ding、Zhiyu Luo、Xing Zhang 的顺序输入作者，指定 Xing Zhang 为唯一通讯作者。共一字段或说明必须按编辑部回复处理。
4. 按页面文件类别上传正文、作者信息页、介绍信、主图、补充表和补充图。如果系统要求标题页与正文合并，使用可编辑来源合并生成，不要遗漏资助和利益冲突。不要上传中文审阅稿、内部说明或脚本。
5. 从官方作者指南下载当前 ACP 表，核对姓名及顺序，全体作者签署；同时上传机构伦理批准/豁免或编辑部认可的替代材料。不得生成作者签名或审批编号。
6. 完成资金、利益冲突、数据可得性及系统其他声明。写作工具相关问题按真实使用情况回答，不填写与实际工作不符的“未使用”。本次核查发现 JCPSP 明确禁止 AI 编辑/润色；当前稿不满足该项政策。不能勾选与实际工作相反的声明来完成提交。
7. 预览系统生成的合稿 PDF，核对上下标、特殊字符、图像、作者顺序和附件完整性；确认摘要数字与正文一致。由通讯作者完成最终确认及提交。
8. 保存稿件编号、确认邮件及合稿 PDF。只有收到系统确认才算提交成功。费用按官网或正式账单核实，不使用旧推荐表中的约数付款。
9. 收到返修意见后逐条回复，并保留修改痕迹版和清洁版。撤稿或转投需先完成当前期刊流程，避免同时在投。

## 五、仍需补齐的个人信息

Tao Li、Jiacai Li、Zhiyu Luo：邮箱和电话；Yajie Ding：电话。Xing Zhang 的电话及邮箱已有。所有作者的学位/专业资格需核实。以投稿表要求决定是否填写 ORCID，不能代为生成。

## 六、关于期刊指标与录用机会

没有将旧表中的影响因子、分区、预警名单或“命中率最高”写入投稿材料。这些信息不构成录用保证。当前重点是解决期刊与三位共一署名的兼容性，以及公开数据研究的提交资格；若不兼容，应保留真实作者贡献并改选期刊。

核查来源：JCPSP 官方作者指南及官方 FAQ（https://www.cpsp.edu.pk/jcpsp.pk/faqs.php）。官网文字可能更新，最终以投稿日要求及编辑部书面回复为准。
'''
supp='# Supplementary figures\n\nOriginal manuscript Figure 1 is Supplementary Figure S1; original Figure 5 is Supplementary Figure S2. Main Figures 1–3 correspond to original Figures 2–4.\n\n'
for n,old,path in [(1,1,'20261024_03_support/20261024_Figure1.png'),(2,5,'20261025_04_review_support/20261025_Figure5.png')]:
 supp+=f'![]({root/path})\n\n**Supplementary Figure S{n}:** '+captions[old]+'\n\n'
texts={'02_Title_page':titlepage,'03_Cover_letter_draft':cover,'06_Supplementary_figures':supp,'07_投稿步骤与待完成事项':guide,'08_Pre_submission_enquiry':inquiry}
pandoc='/home/zhiyu/miniconda3/envs/test/bin/pandoc'
for name,text in texts.items():
 p=src/(name+'.md');p.write_text(text)
 subprocess.run([pandoc,str(p),'--from','markdown+superscript','-o',str(out/(name+'.docx'))],check=True)
 d=Document(out/(name+'.docx'))
 for sec in d.sections:
  sec.page_width=Inches(8.27);sec.page_height=Inches(11.69)
  sec.top_margin=sec.bottom_margin=sec.left_margin=sec.right_margin=Inches(1)
 chinese=name.startswith('07')
 for st in ['Normal','Body Text']:
  stl=d.styles[st];stl.font.name='Times New Roman';stl.font.size=Pt(11 if chinese else 12);stl.paragraph_format.line_spacing=1.25 if chinese else 2
  if chinese:
   stl.element.get_or_add_rPr().rFonts.set(qn('w:eastAsia'),'Noto Sans CJK SC')
 for shape in d.inline_shapes:
  scale=min(1,Inches(6.2)/shape.width,Inches(7.5)/shape.height);shape.width=int(shape.width*scale);shape.height=int(shape.height*scale)
 d.core_properties.author='';d.core_properties.last_modified_by=''
 d.save(out/(name+'.docx'))
shutil.copy2(root/'20261025_22_article_中文版本合并核对.docx',out/'04_中文完整审阅稿_原图编号.docx')
shutil.copy2(root/'20261025_06_supplementary_material.zip',out/'05_Supplementary_tables.zip')
for new,old,path in [(1,2,'20261024_03_support/20261024_Figure2.png'),(2,3,'20261024_03_support/20261024_Figure3.png'),(3,4,'20261025_04_review_support/20261025_Figure4.png')]:
 shutil.copy2(root/path,out/'09_figures'/f'0{new}_Figure_{new}.png')
(out/'00_先读我.txt').write_text('JCPSP 投稿准备包，尚未提交，不建议直接投。官网禁止AI编辑/润色，与实际修改过程冲突。\n请先阅读07投稿步骤；三位共一与期刊规则冲突、伦理豁免文件及其他待办未解决前，不可作为已完成投稿终稿直接提交。\n保留原始数据结果与既定作者顺序。英文主图3幅，原图1和5转入补充图。中文稿沿用原5图编号。\n包含介绍信及编辑部预询问信，两封信均未发送。\n续接任务：先处理共一兼容性和公开数据研究资格，再补实际研究日期、作者资质与联系信息、近期文献及官方签字表。\n')
body=main.split('# Introduction')[1].split('# Patient Consent')[0]
(out/'10_editable_sources'/'00_build_checks.json').write_text(json.dumps({'source':'20261025_21_article_version_reconciled.md','main_body_words_approx':len(body.split()),'references':24,'main_figures':3,'supplementary_figures':2,'author_order':['Tao Li','Jiacai Li','Yajie Ding','Zhiyu Luo','Xing Zhang'],'submitted':False,'email_sent':False},indent=2))
print('Created',out)

# Export the two supplementary figures as standalone files as well.
shutil.copy2(root/'20261024_03_support/20261024_Figure1.png',out/'09_figures/04_Supplementary_Figure_S1.png')
shutil.copy2(root/'20261025_04_review_support/20261025_Figure5.png',out/'09_figures/05_Supplementary_Figure_S2.png')

# Restore the combined manuscript from its authoritative editable source.
subprocess.run([pandoc, str(src/"01_Manuscript_JCPSP_preparation.md"), "--from", "markdown+superscript", "-o", str(out/"01_Manuscript_JCPSP_preparation.docx")], check=True)
