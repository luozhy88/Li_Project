from pathlib import Path
import re

root = Path(__file__).parents[2]
s = (root / 'edit.Paper/output/20261024_24_article_updated.md').read_text()
s = s.replace('Network pharmacology of core herbs selected from Xiao Ying Tang: thyroid-cancer-associated targets and functional enrichment', 'Network pharmacology and independent transcriptome validation of core herbs selected from Xiao Ying Tang in thyroid cancer')
s = s.replace('revision 15.', 'revision 22.')
s = s.replace('This work was supported by the construction project of the inheritance studio of Wang Wenchun, a nationally renowned expert in traditional Chinese medicine, under the National Administration of Traditional Chinese Medicine (Document No. Guo Zhong Yi Yao Ren Jiao Fa [2010] 59); and Shanghai Clinical Research Center of Traditional Chinese Medicine Oncology Science and Technology Commission of Shanghai Municipality (21MC1930500).', 'This work was supported by the Shanghai Clinical Research Center of Traditional Chinese Medicine Oncology, Science and Technology Commission of Shanghai Municipality (21MC1930500).')
s = s.replace(', and Pengfei Yang^6^', '')
s = s.replace("^6^Clinical College of Chinese Medicine, Gansu University of Chinese Medicine, No. 35 Dingxi East Road, Chengguan District, Lanzhou City, Gansu Province, People's Republic of China\n\n", '')
s = s.replace(' Pengfei Yang: Supervision.', '')
s = s.replace('the unvalidated ancillary docking is reported with its limitations.', 'the unvalidated ancillary docking is reported with its limitations. The in-text citation numbers for the GSE60542 source study and limma were corrected and the corresponding reference entry completed. The exploratory docking pose figure (Figure 4) and the transcriptome-validation composite figure (Figure 5) were added to the main text; the bacterial 4P5A panel was excluded. Corresponding authorship was changed to Xing Zhang.')
s = s.replace('Correspondence: Pengfei Yang, Clinical College of Chinese Medicine, Gansu University of Chinese Medicine, No. 35 Dingxi East Road, Chengguan District, Lanzhou City, Gansu Province, People\'s Republic of China. Tel: +8618293195696; Email: ypf@gszy.edu.cn', 'Correspondence: Xing Zhang, Department of Respiratory Disease, Shuguang Hospital Affiliated to Shanghai University of Traditional Chinese Medicine, Shanghai 201203, China. Tel: 021-20256375; Email: zhangxing2012@shutcm.edu.cn')
s = s.replace('Additional author contact details: Yajie Ding, Email: dingyajie2018@163.com; Xing Zhang, Email: zhangxing2012@shutcm.edu.cn, Tel: 021-20256375.', 'Additional author contact details: Yajie Ding, Email: dingyajie2018@163.com.')

# Add the two references at first-use order and shift the original references 8–20 once.
# This renumbering must run BEFORE inserting new text that already uses the final
# citation numbers (^8^, ^9^); otherwise the inserted citations are shifted again.
s = re.sub(r'\^([0-9,]+)\^', lambda m: '^' + ','.join(str(int(x)+2 if int(x) >= 8 else int(x)) for x in m.group(1).split(',')) + '^', s)
s = re.sub(r'(?m)^(\d+)\. ', lambda m: f'{int(m.group(1))+2}. ' if int(m.group(1)) >= 8 else m.group(0), s)
s = s.replace('10. Yu G, Wang LG, Han Y, He QY.', '8. Tarabichi M, Saiselet M, Trésallet C, Hoang C, Larsimont D, Andry G, et al. Revisiting the transcriptional analysis of primary tumours and associated nodal metastases with enhanced biological and statistical controls: application to thyroid cancer. Br J Cancer. 2015;112(10):1665–74. doi:10.1038/bjc.2014.665.\n\n9. Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. 2015;43(7):e47. doi:10.1093/nar/gkv007.\n\n10. Yu G, Wang LG, Han Y, He QY.')

a = s.index('# Abstract'); b = s.index('# Introduction')
s = s[:a] + '''# Abstract

**Objective:** To explore thyroid-cancer-associated targets and pathways of core herbs selected from Xiao Ying Tang (XYT), with Zhe Bei Mu (ZBM) and Xia Ku Cao (XKC) as individual-herb comparators, and to assess these prespecified gene sets in an independent thyroid-cancer transcriptome.

**Study design:** Descriptive in silico network pharmacology study with independent paired transcriptome validation and ancillary docking assessment.

**Methodology:** Database-derived compound–gene associations were intersected with a GeneCards-derived disease set. GO and KEGG over-representation analyses used Benjamini–Hochberg adjustment. GSE60542/GPL570 was analysed using the submitter-normalized fRMA matrix. Primary papillary thyroid cancer and normal thyroid samples were restricted to 24 patients with both tissue types; repeated blocks were averaged within patient and tissue. A paired limma model was fitted across all measured genes, followed by fixed-background Fisher tests for the prespecified XYT, XKC, and ZBM sets.

**Results:** The disease-associated intersections contained 78, 842, and 1,595 genes for ZBM, XKC, and XYT, respectively, with the individual-herb intersections nested within XYT. The transcriptome analysis measured 20,848 genes and identified 666 genes at FDR <0.05 and |log2FC| ≥1. XYT was enriched among differentially expressed genes in the all-measured background (odds ratio 2.67, FDR =6.08 × 10^-16^) and disease-measured background (1.80, 5.31 × 10^-7^). XKC showed corresponding enrichment (3.05, 2.44 × 10^-13^; 2.08, 5.31 × 10^-7^), whereas ZBM was not significant. POLB, CACNA1H, and CYP3A4 did not meet the prespecified differential-expression threshold.

**Conclusion:** The independent transcriptome supports differential-expression enrichment of the prespecified XYT and XKC sets, but not ZBM. These results provide transcript-level association evidence rather than proof of drug action, target engagement, pathway direction, or clinical efficacy.

**Keywords:** Thyroid Neoplasms; Medicine, Chinese Traditional; Gene Ontology; Computational Biology; Transcriptome

''' + s[b:]

anchor = 'Gene symbols were mapped to Entrez identifiers with org.Hs.eg.db'
method = '''Independent transcriptome validation used the prespecified GEO series GSE60542 on GPL570, whose original study examined papillary thyroid carcinoma and normal thyroid and other tissue types.^8^ The downloaded series matrix was the submitter-provided fRMA-normalized matrix; values were treated as log2 expression and were not re-normalized. Only primary papillary thyroid cancer and normal thyroid samples were retained. Repeated blocks from the same patient and tissue type were averaged, giving 24 patients and 48 patient-by-tissue observations. GPL570 probes with one non-empty gene symbol were retained; multiple probes per gene were averaged on the log2 scale. Differential expression used limma with a patient fixed-effect design (`~ patient_id + tissue_group`), normal thyroid as reference, and positive log2 fold change indicating higher tumour expression.^9^ The FDR was calculated across all 20,848 measured genes; the threshold was FDR <0.05 and absolute log2 fold change ≥1. Enrichment of the fixed XYT, XKC, and ZBM sets used one-sided Fisher tests against all measured genes and the archived disease set intersected with measured genes. The six tests were adjusted together by Benjamini–Hochberg. This analysis did not redefine candidate sets or thresholds.

'''
s = s.replace(anchor, method + anchor)

anchor = 'β-Sitosterol (CID 222284) was represented in all three ingredient sets'
result = '''**Independent transcriptome validation.** The final paired analysis included 24 patients, 28 primary-tumour blocks, and 27 normal-thyroid blocks, aggregated to 48 observations. The model had design rank 25 and 23 residual degrees of freedom. Of 20,848 measured genes, 666 met FDR <0.05 and |log2FC| ≥1; 347 were higher and 319 lower in tumour than normal thyroid. XYT contained 111 differentially expressed genes among 1,518 measured genes and was enriched against all measured genes (odds ratio 2.67, P =1.01 × 10^-16^, FDR =6.08 × 10^-16^) and the disease-measured background (1.80, P =3.16 × 10^-7^, FDR =5.31 × 10^-7^). XKC contained 69 differentially expressed genes among 805 measured genes and was enriched in both backgrounds (odds ratios 3.05 and 2.08; FDR =2.44 × 10^-13^ and 5.31 × 10^-7^). ZBM contained 4 differentially expressed genes among 77 measured genes; enrichment was not significant (FDR =0.278 and 0.482). POLB, CACNA1H, and CYP3A4 were measured but did not meet the differential-expression threshold (log2FC =−0.070, −0.053, and −0.031, respectively). These results do not establish pharmacological direction or direct binding.

'''
s = s.replace(anchor, result + anchor)

anchor = 'Experimental follow-up should use chemically characterised preparations'
discussion = '''The independent GSE60542 analysis provides a transcript-level check of the prespecified database-derived sets. XYT and XKC were enriched among genes differing between primary papillary thyroid cancer and normal thyroid, while ZBM was not. This is not a ranking of herb efficacy: the sets are nested, the candidates were defined from database and disease-gene resources, and the transcriptome comparison is observational. Expression differences do not demonstrate target engagement, reversal of tumour state, or patient benefit. Replication in independent cohorts and functional experiments are required.

'''
s = s.replace(anchor, discussion + anchor)
s = s.replace('The study did not measure constituent concentrations, bioavailability, target engagement, cellular responses, or clinical outcomes.', 'The transcriptome validation used one public microarray cohort, a modest paired sample, and the submitter’s fRMA processing; it is not a treatment comparison and cannot establish causality or therapeutic response. The study did not measure constituent concentrations, bioavailability, target engagement, cellular responses, or clinical outcomes.')

a = s.index('# Conclusion'); b = s.index('# Conflict of Interest')
s = s[:a] + '''# Conclusion

The selected XYT core herbs yielded a larger disease-associated gene set that contained the individual-herb intersections. Functional annotation highlighted oxidative stress, calcium regulation, and signalling. In an independent paired GSE60542 transcriptome, the prespecified XYT and XKC sets were enriched among differentially expressed genes, whereas ZBM was not. POLB, CACNA1H, and CYP3A4 did not meet the prespecified gene-level threshold. These findings provide transcript-level association evidence for prioritising experiments, but do not establish direct binding, pharmacological direction, synergy, or clinical efficacy.

''' + s[b:]
s = s.replace('The reported work analysed database records and molecular structures without recruiting human participants or using animals.', 'The reported work analysed database records, molecular structures, and a public de-identified GEO transcriptome without new recruitment, intervention, or identifiable patient data. No human participants or animals were newly studied.')
s = s.replace('Analysis results, gene lists, and supporting tables are provided as supplementary material. The underlying data are available from the corresponding author upon reasonable request. The analysis code is not publicly available.', 'Analysis results, gene lists, sample and pairing manifests, full differential-expression results, enrichment tests, and QC figures are provided as supplementary material. The transcriptome source data are publicly available from GEO under accession GSE60542; the underlying network-pharmacology source data are available from the corresponding author upon reasonable request. The analysis code is retained locally for audit and is not included in the public supplementary package.')

s = s.replace('−9.68 kcal/mol for CaV3.2 (9AYG). The 4P5A output was excluded from human-target interpretation.', '−9.68 kcal/mol for CaV3.2 (9AYG) (Figure 4). The 4P5A output was excluded from human-target interpretation.')
anchor = 'No reference-ligand redocking assessment or direct target-binding experiment was available.'
fig = '''No reference-ligand redocking assessment or direct target-binding experiment was available.

![](20261024_03_support/20261024_Figure4.png)

**Figure 4:** Exploratory docking poses of β-sitosterol (CID 222284). (A) Human CYP3A4 (PDB 6UNI); (B) human CaV3.2/CACNA1H (PDB 9AYG). Left panels show the local view of the docked ligand (green sticks) and surrounding residues; right panels show the global receptor view. Poses correspond to the minimum supplementary smina scores of −9.12 and −9.68 kcal/mol, respectively. The bacterial 4P5A entry was excluded from human-target interpretation (see Methodology). These poses are exploratory structural outputs and are not validated binding evidence.

'''
s = s.replace(anchor, fig)

s = s.replace('Zhiyu Luo is an employee of Jiangsu Sanshu Biotechnology Co., Ltd. The authors declare no conflicts of interest related to this work.', 'Zhiyu Luo is an employee of Jiangsu Sanshu Biotechnology Co., Ltd. The remaining authors declare no conflicts of interest related to this work.')

s = s.replace('347 were higher and 319 lower in tumour than normal thyroid.', '347 were higher and 319 lower in tumour than normal thyroid (Figure 5).')
anchor = 'These results do not establish pharmacological direction or direct binding.'
fig5 = '''These results do not establish pharmacological direction or direct binding.

![](20261024_03_support/20261024_Figure5.png)

**Figure 5:** Independent transcriptome validation in GSE60542 (24 paired primary papillary thyroid carcinomas vs matched normal thyroid). (A) Volcano plot of the paired limma analysis across 20,848 measured genes; dashed lines mark FDR =0.05 and |log2FC| =1. Coloured points denote measured genes in the prespecified ZBM (orange), XKC (teal), and XYT (purple) sets. (B) Odds ratios with two-sided 95% confidence intervals for enrichment of each prespecified set among differentially expressed genes, against the all-measured and disease-measured backgrounds; significance was assessed with the prespecified one-sided Fisher tests, and FDR values are from Benjamini–Hochberg adjustment across the six tests.

'''
s = s.replace(anchor, fig5)

out = root / 'edit.Paper/output/20261025_07_article_funding_updated.md'
out.write_text(s)
print(out)
