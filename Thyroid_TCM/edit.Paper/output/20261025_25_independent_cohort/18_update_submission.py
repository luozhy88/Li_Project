from pathlib import Path
import re,shutil,subprocess,zipfile,json,hashlib
from docx import Document
from docx.shared import Inches,Pt
R=Path('/home/data1/zhiyu/data/Github/Li_Project/Thyroid_TCM/edit.Paper/output');A=R/'20261025_25_independent_cohort';P=R/'20261025_19_JCPSP_submission';S=R/'20261025_24_投稿整理归档/10_editable_sources'
backup=A/'17_before_update';backup.mkdir(exist_ok=True)
for n in ['01_Manuscript_JCPSP_preparation.docx','03_Cover_letter_draft.docx','05_Supplementary_tables.zip','06_Supplementary_figures.docx']:
 if not (backup/n).exists():shutil.copy2(P/n,backup/n)
s=(A/'17_before_update/01_manuscript_source.md').read_text()
def replace_paragraph(prefix,new):
 global s
 old=[p for p in s.split('\n\n') if p.startswith(prefix)];assert len(old)==1,prefix;s=s.replace(old[0],new,1)
abstract='''# Abstract

**Objective:** To characterise genes associated with core herbs selected from Xiao Ying Tang (XYT) and assess their expression across two independent papillary thyroid carcinoma (PTC) cohorts.

**Study design:** Secondary analysis of public databases and paired tumour–normal expression data.

**Methodology:** XYT-, Zhe Bei Mu (ZBM)- and Xia Ku Cao (XKC)-associated genes were intersected with a GeneCards disease list and analysed for functional enrichment. Paired limma models assessed 24 patients in GSE60542 and 44 in GSE33630. Differential expression required Benjamini–Hochberg-adjusted P <0.05 and |log2 fold change| ≥1. Candidate-set enrichment was tested against all measured genes and disease-associated measured genes.

**Results:** Disease intersections contained 78 ZBM, 842 XKC and 1,595 XYT genes. Among 20,848 measured genes per cohort, 666 and 745 met the differential-expression criteria. XYT and XKC were enriched in both cohorts. Against disease-restricted backgrounds, odds ratios were 1.80 and 2.08 in GSE60542, and 1.59 and 1.80 in GSE33630; all adjusted P values were <0.001. Of the 666 initial differential-expression findings, 664 had concordant directions and 536 met the same criteria in GSE33630. ZBM enrichment was not significant. POLB, CACNA1H and CYP3A4 did not meet the combined differential-expression criteria in either cohort.

**Conclusion:** XYT- and XKC-associated genes showed reproducible enrichment among genes altered in PTC tissue. These disease associations prioritise candidates for functional investigation.

**Keywords:** Thyroid Neoplasms; Medicine, Chinese Traditional; Gene Ontology; Computational Biology; Transcriptome

'''
a=s.index('# Abstract');b=s.index('# Introduction');s=s[:a]+abstract+s[b:]
s=s.replace('then assessed the fixed candidate sets in paired papillary thyroid carcinoma (PTC) and normal thyroid samples.','then assessed the fixed candidate sets in paired papillary thyroid carcinoma (PTC) and normal thyroid samples from two independent cohorts.')
replace_paragraph('The external expression analysis used', '''The expression analyses used two Affymetrix GPL570 cohorts: GSE60542 from Belgium and France,^10^ and GSE33630 from the Ukrainian Chernobyl Tissue Bank.^NEW^ Submitted fRMA- and RMA-normalised log2 matrices, respectively, were analysed separately. Only matched primary PTC and non-tumoral thyroid samples were included. Repeated GSE60542 tissue blocks were averaged within patient and tissue type. GSE33630 title-based pairing was confirmed against clinical bank identifiers; five unpaired tumours, one unpaired normal sample and 11 anaplastic cancers were excluded. The cohorts shared no GEO sample identifiers. Probes with missing or ambiguous gene annotations were excluded; probes mapping to one gene were averaged on the log2 scale, yielding 20,848 genes per cohort. Paired limma models (`~ patient_id + tissue_group`) used R 4.3.1, limma 3.58.1 and robust empirical Bayes moderation.^11^ Positive log2 fold changes indicate higher tumour expression. Differential expression required Benjamini–Hochberg-adjusted P <0.05 and |log2FC| ≥1. The unchanged XYT, XKC and ZBM sets were tested using one-sided Fisher tests against all measured genes and the 10,655 disease-associated measured genes. Six P values were adjusted within each cohort. Conditional odds ratios have two-sided exact 95% confidence intervals. Cross-cohort agreement was assessed by Spearman correlation, expression direction and repeat fulfilment of the differential-expression criteria. Exposure-stratified GSE33630 sensitivity analyses used the same paired model and thresholds; 12 enrichment P values were adjusted together.''')
# Add the new results before the exploratory docking paragraph.
newresults='''GSE33630 contributed 44 matched pairs and identified 745 differentially expressed genes (407 higher and 338 lower in tumour tissue). XYT included 112/1,518 measured candidates and XKC 68/805. Their odds ratios were 2.35 and 2.64 against all measured genes, and 1.59 (95% CI 1.27–1.97) and 1.80 (1.36–2.35) against disease-associated genes (adjusted P =4.99 × 10^-5^ for both disease-restricted tests). ZBM included 2/77 genes and was not enriched (adjusted P =0.913 for both backgrounds). Across 20,848 genes, log2 fold changes correlated between cohorts (Spearman ρ =0.786). Of 666 initial differential-expression findings, 664 agreed in direction and 536 met the same criteria in GSE33630, including 83/111 XYT and 51/69 XKC genes (Supplementary Figure S3). POLB, CACNA1H and CYP3A4 again failed the combined criteria. XYT and XKC enrichment persisted in both exposure strata (22 pairs each; all adjusted P <0.001), whereas ZBM enrichment remained non-significant.'''
s=s.replace('β-Sitosterol (CID 222284) was common',newresults+'\n\nβ-Sitosterol (CID 222284) was common',1)
s=s.replace('In GSE60542, XYT and XKC candidates were over-represented among genes altered in tumours,','In both cohorts, XYT and XKC candidates were over-represented among genes altered in tumours,',1)
s=s.replace('The expression cohort represents PTC, and differences in tissue composition may contribute to the observed changes.','Both cohorts represent PTC, and tissue composition may contribute to expression differences. GSE33630 comprised younger patients (median age 16.82 years; range 12.38–24.55), with equal numbers classified as exposed and unexposed to radiation. Agreement across these cohorts supports reproducibility of the disease associations; it does not demonstrate a response to herbal treatment.',1)
s=s.replace('XYT and XKC candidate genes were over-represented among genes altered in paired PTC tissue, including against a disease-restricted background.','XYT and XKC candidate genes were over-represented among genes altered in paired PTC tissue in two independent cohorts, including against disease-restricted backgrounds.',1)
s=s.replace('Expression data and platform annotations are available from GEO under GSE60542 and GPL570.','Expression data and platform annotations are available from GEO under GSE60542, GSE33630 and GPL570. Cross-cohort comparisons and exposure-stratified results are included in the supplementary material.')
# Insert the new source reference in first-citation order.
body,refs=s.split('# References\n',1)
body=re.sub(r'\^(\d+(?:,\d+)*)\^',lambda m:'^'+','.join(str(int(x)+1 if int(x)>=11 else int(x)) for x in m[1].split(','))+'^',body)
body=body.replace('^NEW^','^11^')
refs=re.sub(r'^(\d+)\.',lambda m:str(int(m[1])+1 if int(m[1])>=11 else int(m[1]))+'.',refs,flags=re.M)
newref='11. Dom G, Tarabichi M, Unger K, Thomas G, Oczko-Wojciechowska M, Bogdanova T, et al. A gene expression signature distinguishes normal tissues of sporadic and radiation-induced papillary thyroid carcinomas. Br J Cancer. 2012;107(6):994–1000. doi:10.1038/bjc.2012.302.\n\n'
refs=refs.replace('12. Ritchie',newref+'12. Ritchie',1);s=body+'# References\n'+refs
(A/'19_manuscript_two_cohorts.md').write_text(s)
pandoc='/home/zhiyu/miniconda3/envs/test/bin/pandoc'
def convert(textfile,out):
 subprocess.run([pandoc,str(textfile),'--from','markdown+superscript','-o',str(out)],check=True)
 d=Document(out)
 for sec in d.sections:
  sec.page_width=Inches(8.27);sec.page_height=Inches(11.69);sec.top_margin=sec.bottom_margin=sec.left_margin=sec.right_margin=Inches(1)
 for style in ['Normal','Body Text']:
  d.styles[style].font.name='Times New Roman';d.styles[style].font.size=Pt(12);d.styles[style].paragraph_format.line_spacing=2
 for para in d.paragraphs:
  if para.text=='Abstract':para.paragraph_format.page_break_before=True
 for shape in d.inline_shapes:
  scale=min(1, Inches(6.2)/shape.width,Inches(7.3)/shape.height);shape.width=int(shape.width*scale);shape.height=int(shape.height*scale)
 d.save(out)
convert(A/'19_manuscript_two_cohorts.md',P/'01_Manuscript_JCPSP_preparation.docx')
# Add only the new supplementary figure; preserve the existing S1/S2 figures.
shutil.copy2(A/'16_Supplementary_Figure_S3.png',P/'09_figures/06_Supplementary_Figure_S3.png')
caption='Supplementary Figure S3: Independent replication in GSE33630 (44 matched patients). (A) Paired-limma differential expression: BH-adjusted P <0.05 and |log2FC| ≥1. (B) Candidate-set enrichment against all 20,848 measured genes or 10,655 disease-associated measured genes. Points indicate conditional odds ratios; bars show two-sided exact 95% confidence intervals. q labels are BH-adjusted P values from six one-sided Fisher tests in GSE33630. (C) Gene-level log2 fold changes across cohorts; blue points mark the 666 genes meeting the differential-expression criteria in GSE60542. Of these, 664 had the same direction and 536 met the same criteria with the same direction in GSE33630. The correlation uses all 20,848 genes.'
d=Document(backup/'06_Supplementary_figures.docx');d.add_page_break();d.add_picture(str(A/'16_Supplementary_Figure_S3.png'),width=Inches(6.2));d.add_paragraph(caption);d.save(P/'06_Supplementary_figures.docx')
(A/'20_S3_caption.txt').write_text(caption)
# Update the existing cover letter in place; no individual employment mention.
d=Document(backup/'03_Cover_letter_draft.docx')
for para in d.paragraphs:
 t=para.text
 t=t.replace('an independent papillary thyroid carcinoma expression dataset','two independent papillary thyroid carcinoma expression datasets').replace('from 24 patients, with repeated tissue blocks averaged within each patient and tissue type','from 24 patients in GSE60542 and 44 in GSE33630, analysed separately with paired models').replace('were over-represented among differentially expressed genes under both backgrounds','were over-represented among differentially expressed genes under both backgrounds in both cohorts')
 t=t.replace('The smaller Zhe Bei Mu set was not significantly enriched,','Of the 666 differential-expression findings in GSE60542, 536 met the same criteria with concordant direction in GSE33630. The smaller Zhe Bei Mu set was not significantly enriched,')
 t=t.replace('are disclosed on the title page','are disclosed in the manuscript')
 if t!=para.text:para.text=t
assert 'Zhiyu Luo' not in '\n'.join(x.text for x in d.paragraphs)
d.save(P/'03_Cover_letter_draft.docx')
# Preserve existing supplemental results and add reproducible numerical evidence, without code.
files=[('GSE33630_02_sample_manifest_all.csv','25_GSE33630_sample_manifest.csv'),('11_verified_paired_manifest.csv','26_GSE33630_verified_pairs.csv'),('GSE33630_03_all_measured_gene_DE.csv','27_GSE33630_all_gene_DE.csv'),('GSE33630_enrichment_with_CI.csv','28_GSE33630_set_enrichment.csv'),('GSE33630_cross_cohort_all_genes.csv','29_cross_cohort_gene_comparison.csv'),('GSE33630_concordance_summary.csv','30_cross_cohort_summary.csv'),('GSE33630_exposure_sensitivity.csv','31_exposure_sensitivity.csv'),('GSE33630_Ex_all_gene_DE.csv','32_exposed_stratum_DE.csv'),('GSE33630_NE_all_gene_DE.csv','33_unexposed_stratum_DE.csv'),('GSE33630_07_QC_PCA.pdf','34_GSE33630_PCA.pdf'),('GSE33630_09_session_info.txt','35_GSE33630_session_info.txt'),('GSE33630_11_probe_mapping_qc.csv','36_GSE33630_probe_QC.csv'),('GSE33630_18_model_qc.csv','37_GSE33630_model_QC.csv'),('GSE33630_methods.txt','38_GSE33630_methods.txt')]
with zipfile.ZipFile(backup/'05_Supplementary_tables.zip') as zin,zipfile.ZipFile(P/'05_Supplementary_tables.zip','w',zipfile.ZIP_DEFLATED) as zout:
 for name in zin.namelist():
  data=zin.read(name)
  if name=='00_README.txt':data+=b'\n\nAdditional independent cohort GSE33630: files 25-38. Primary analysis uses 44 matched patients; 11 ATC and 6 unpaired tissues are excluded. Files 29-30 compare gene effects with GSE60542; file 31 reports exposure sensitivity, BH across 12 tests. Dataset-specific gene and Fisher P values are adjusted separately. Code is not included.\n'
  zout.writestr(name,data)
 for source,dest in files:zout.write(A/source,dest)
# Verification against actual Word contents and the original submission author metadata.
t='\n'.join(x.text for x in Document(P/'01_Manuscript_JCPSP_preparation.docx').paragraphs)
assert all(x in t for x in ['Xing Zhang','zhangxing2012@shutcm.edu.cn','GSE33630','745','536','21MC1930500','ZX0551-25','2025JZ02'])
assert not any(x in t for x in ['Pengfei Yang','Wang Wenchun','2024_09_05'])
assert len(Document(P/'01_Manuscript_JCPSP_preparation.docx').inline_shapes)==3
assert len(Document(P/'06_Supplementary_figures.docx').inline_shapes)==3
assert len(re.findall(r'^\d+\.',refs,re.M))==25
bodycount=re.sub(r'!\[\].*|\*\*Figure [^\n]+','',s.split('# Introduction')[1].split('# Conflict of Interest')[0])
checks={'abstract_words':len(abstract.split('**Keywords:**')[0].split())-2,'body_words_approx':len(bodycount.split()),'references':25,'main_figures':3,'supplementary_figures':3,'authors':5,'main_sha256':hashlib.sha256((P/'01_Manuscript_JCPSP_preparation.docx').read_bytes()).hexdigest()}
(A/'21_manuscript_checks.json').write_text(json.dumps(checks,indent=2));print(checks)
shutil.make_archive(str(P),'zip',P.parent,P.name)
