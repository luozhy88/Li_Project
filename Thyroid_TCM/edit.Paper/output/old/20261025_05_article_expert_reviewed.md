# Core herbs selected from Xiao Ying Tang: network pharmacology and paired transcriptome analysis in papillary thyroid cancer

Tao Li^1^\*, Jiacai Li^2^\*, Yajie Ding^3^\*, Zhiyu Luo^4^, Xing Zhang^5^, and Pengfei Yang^6^

^1^Department of Acupuncture and Moxibustion, Affiliated Hospital of Hangzhou Normal University, Hangzhou, Zhejiang, China.

^2^School of Clinical Medicine, Shanghai University of Medicine & Health Sciences, Shanghai, China.

^3^Department of Oncology, Longhua Hospital Affiliated to Shanghai University of Traditional Chinese Medicine, Shanghai, China.

^4^Jiangsu Sanshu Biotechnology Co., Ltd., No. 259 Shibei Gaoxin Road, Tangzha Town Subdistrict, Chongchuan District, Nantong, Jiangsu, China.

^5^Department of Respiratory Disease, Shuguang Hospital Affiliated to Shanghai University of Traditional Chinese Medicine, Shanghai 201203, China.

^6^Clinical College of Chinese Medicine, Gansu University of Chinese Medicine, No. 35 Dingxi East Road, Chengguan District, Lanzhou City, Gansu Province, People's Republic of China

\*Tao Li, Jiacai Li, and Yajie Ding contributed equally to this work and share first authorship.

Correspondence: Pengfei Yang, Clinical College of Chinese Medicine, Gansu University of Chinese Medicine, No. 35 Dingxi East Road, Chengguan District, Lanzhou City, Gansu Province, People's Republic of China. Tel: +8618293195696; Email: ypf@gszy.edu.cn

Additional author contact details: Yajie Ding, Email: dingyajie2018@163.com; Xing Zhang, Email: zhangxing2012@shutcm.edu.cn, Tel: 021-20256375.

# Abstract

**Objective:** To characterise thyroid-cancer-associated genes linked to core herbs selected from Xiao Ying Tang (XYT) and assess their expression in an external papillary thyroid carcinoma (PTC) cohort.

**Methods:** Compound-associated genes for XYT, Zhe Bei Mu (ZBM), and Xia Ku Cao (XKC) were intersected with a GeneCards disease list and analysed for Gene Ontology and KEGG enrichment. In GSE60542, repeated tissue blocks were averaged within patient and tissue type. Tumour–normal differences were analysed with paired limma models for 24 patients. Differential expression required a Benjamini–Hochberg-adjusted P value <0.05 and an absolute log2 fold change ≥1. Candidate-set enrichment was tested against all measured genes and a background restricted to disease-associated measured genes.

**Results:** The disease intersections comprised 78 ZBM, 842 XKC, and 1,595 XYT genes; both individual-herb sets were contained within XYT. Oxidative-stress and calcium-related annotations were enriched. Of 20,848 genes assessed in GSE60542, 666 met the differential-expression criteria. These included 111/1,518 measured XYT genes, 69/805 XKC genes, and 4/77 ZBM genes. XYT and XKC remained enriched against the disease-restricted background (odds ratios 1.80 and 2.08; adjusted P =5.31 × 10^-7^ for each). ZBM enrichment was not statistically significant. None of POLB, CACNA1H, or CYP3A4 met the combined differential-expression criteria.

**Conclusion:** Genes associated with XYT and XKC were over-represented among genes altered in PTC tissue. This supports their relevance for further biological investigation, but does not establish a treatment effect or validate the proposed compound–protein interactions.

**Keywords:** Thyroid Neoplasms; Medicine, Chinese Traditional; Gene Ontology; Computational Biology; Transcriptome

# Introduction

Thyroid cancer includes malignancies with distinct molecular features and clinical courses.^1^ Surgery, radioactive iodine in selected patients, and thyroid hormone management form the basis of differentiated thyroid cancer care.^2^ Persistent or treatment-resistant disease has prompted investigation of additional therapeutic approaches. Because molecular drivers differ across thyroid cancer subtypes, candidate mechanisms identified from broad disease databases need to be examined in a defined tumour context.^3^

Network pharmacology can link herbal constituents to genes and biological processes. Its value lies in generating testable hypotheses, although database coverage and predicted interactions can strongly influence the resulting gene sets. An external tumour-expression dataset can assess whether those genes are altered in disease tissue, while leaving their response to herbal treatment unresolved.

We studied core herbs selected from Xiao Ying Tang for their relevance to thyroid disease; XYT refers to this selected set throughout. Zhe Bei Mu (ZBM; Fritillariae Thunbergii Bulbus) and Xia Ku Cao (XKC; Prunellae Spica) were examined separately to describe their contributions to the broader annotation profile. We compared compound-associated gene coverage and functional annotations, then assessed the fixed candidate sets in paired papillary thyroid carcinoma (PTC) and normal thyroid samples. Archived docking results for β-sitosterol were also reviewed as exploratory structural observations.

# Methodology

This computational study was conducted at Shanghai University of Medicine & Health Sciences, Shanghai, China. It combined archived herb–compound and compound–gene records, functional enrichment, a secondary analysis of public tumour-expression data, and an assessment of archived docking outputs. No participants were recruited and no treatment experiments were performed. Some retrieval dates and software versions from the original network-pharmacology analysis were not recoverable; the newly performed transcriptome analysis is documented separately below.

Herb–ingredient records from TCMSP were matched to PubChem compound identifiers (CIDs) through lower-case compound names and a synonym table.^4,5^ The integrated compound–gene table included BindingDB, BATMAN-TCM known and predicted interactions, IUPHAR/BPS Guide to PHARMACOLOGY, and TCMSP target records.^6,7,8^ The build script referenced a September 2024 BindingDB file and a UniProt mapping file dated 5 September 2024. BindingDB, BATMAN-TCM known, and Guide to PHARMACOLOGY records were assigned an inclusion value of 100. BATMAN-TCM prediction scores were multiplied by 100; the implemented `Gene_Value > 80` filter therefore retained predictions with an original score >0.80. These assigned values do not represent a common experimental-confidence scale. TCMSP target records had a differently capitalised score field (`Gene_value`) and missing `Gene_Value` values after merging, so their target records were excluded by the filter. TCMSP remained the source of herb–ingredient records. No oral-bioavailability or drug-likeness threshold was applied. BindingDB records were not filtered by binding affinity or organism before gene-symbol integration.

The XYT input contained the herb labels E Zhu, Bai Shao, Dang Gui, Chai Hu, Qing Pi, Huang Yao Zi, Xia Ku Cao, Zao Jiao Ci, Mu Dan Pi, Bai Zhu, Fu Ling, Bai Mao Xia Ku Cao, Xuan Shen, and Zhe Bei Mu. This inventory defines the selected set analysed here, rather than the complete prescription. ZBM and XKC records were extracted using their Chinese herb names and processed by the same workflow. Gene symbols were deduplicated within each CID, and distinct CIDs were counted as mapped ingredients, ingredients linked to retained genes, and ingredients linked to disease-overlapping genes. Duplicate herb–compound records were collapsed within each set. Herb doses and constituent concentrations were unavailable.

The GeneCards export contained 12,987 unique gene symbols.^9^ Its filename, `GeneCards-SearchResults_Thyroid_tumor.csv`, indicates the search topic, but the exact historical query syntax and retrieval date could not be recovered. The investigators confirmed that no filtering was applied before export. The full list was used without a relevance-score threshold. Ingredient-associated gene symbols were split on semicolons, deduplicated, and intersected with this list. Reconstructing the ingredient-associated sets from source-labelled records reproduced the archived sets exactly.

The external expression analysis used GSE60542, an Affymetrix GPL570 study of PTC and associated tissues.^10^ We used the submitter's fRMA-normalised log2 expression matrix without further transformation or normalisation. Only primary PTC and normal thyroid samples were eligible; lymph-node tissue and nodal metastases were excluded. Patients with both eligible tissue types were retained. Repeated blocks were averaged within each patient and tissue type. Probes with missing, placeholder, or ambiguous multiple-gene annotations were excluded; retained probes mapping to the same gene were averaged on the log2 scale. This yielded 20,848 genes from 42,904 of 54,675 probes. Differential expression was analysed in R 4.3.1 with limma 3.58.1 using `~ patient_id + tissue_group`, normal thyroid as the reference, and robust empirical Bayes moderation.^11^ Positive log2 fold changes indicate higher tumour expression. P values were adjusted across all measured genes by the Benjamini–Hochberg method; differential expression required adjusted P <0.05 and |log2FC| ≥1. The XYT, XKC, and ZBM candidate sets and thresholds were fixed before this analysis. One-sided Fisher exact tests assessed over-representation among differentially expressed genes, using either all 20,848 measured genes or the 10,655 measured genes in the archived disease list as the background. The six P values were adjusted together. Forest plots show conditional odds-ratio estimates and two-sided exact 95% confidence intervals; their P values come from the specified one-sided tests.

For the original functional analysis, unique gene symbols were mapped to Entrez identifiers with `bitr` and org.Hs.eg.db, and all returned identifiers were passed to clusterProfiler.^12^ No custom rule resolved one-to-many mappings. GO analysis included biological process (BP), molecular function (MF), and cellular component (CC); KEGG analysis used human annotations. Both analyses specified Benjamini–Hochberg adjustment, `pvalueCutoff = 0.05`, and `qvalueCutoff = 0.2`. No custom background was supplied. Exported results showed annotation backgrounds of 18,888 genes for BP, 18,522 for MF, 19,894 for CC, and 8,859 for KEGG. Terms with adjusted P <0.05 were reported. Raw hypergeometric P values were checked against the archived counts; the original adjusted P values were retained because the complete pre-filtering test lists were unavailable. Historical package and annotation versions were not archived.

GO figures display the union of the three lowest-adjusted-P terms per set within each ontology. KEGG figures display the union of the five lowest-adjusted-P pathways per set, with the thyroid cancer pathway (hsa05216) added because of the study question. GeneRatio uses the annotated input size recorded for each analysis, which may differ across ontologies. Missing bars indicate that a term was absent from an exported significant-result table. The figures describe annotation profiles; no formal between-herb comparison of enrichment or efficacy was performed.

Ligands were prepared with Open Babel and receptor structures obtained from RCSB PDB.^13,14^ Initial screening used AutoDock Vina; supplementary β-sitosterol docking used smina.^15,16^ The reviewed human structures were CYP3A4 (6UNI) and apo CaV3.2/CACNA1H (9AYG).^14,17^ The bacterial 4P5A entry was excluded from interpretation as human POLB. The initial Vina script specified a 20 × 20 × 20 Å box centred on the mean receptor ATOM coordinates rounded to whole angstroms, exhaustiveness 32, and five CPU threads. Supplementary reports recorded exhaustiveness 8, nine output modes, and CPU-only execution. The corresponding pipeline used one CPU thread and a reference-based automatic box with `--autobox_add 4`; logs confirmed that branch was used. A random seed and scoring-function override were not specified. Historical executable versions and complete engine invocation logs were unavailable. Receptor identity, retained cofactors, and box-reference contents were reviewed before interpreting the archived outputs.

# Results

The ingredient tables contained 16, 68, and 803 distinct non-missing PubChem CIDs for ZBM, XKC, and XYT, respectively. Of these, 8, 57, and 501 had retained gene associations; 8, 57, and 492 were associated with at least one gene in the disease intersection. The corresponding unique ingredient-associated gene totals were 107, 1,094, and 2,131. Of these, 78/107 (72.9%), 842/1,094 (77.0%), and 1,595/2,131 (74.8%) overlapped with the common disease set. Disease-only counts were 12,909, 12,145, and 11,392, while ingredient-only counts were 29, 252, and 536 (Figure 1). All 78 ZBM and all 842 XKC disease-overlapping genes were contained in the 1,595-gene XYT intersection. XYT covered the largest gene set, while the proportions overlapping the disease list were similar across the three sets.

![](20261024_03_support/20261024_Figure1.png)

**Figure 1:** Gene-set intersections for ZBM, XKC, and the XYT core-herb set. Shared disease genes total 78, 842, and 1,595, respectively. All comparisons use the same 12,987-gene disease set. Circle areas are schematic and are not proportional to set sizes. XYT denotes the core herbs selected from Xiao Ying Tang for thyroid disease.

The saved GO result tables contained 1,026, 3,925, and 4,782 terms with adjusted P <0.05 for ZBM, XKC, and XYT. Their annotated input sizes were 78/78/78 for ZBM, 835/840/839 for XKC, and 1,582/1,588/1,587 for XYT in BP/MF/CC order. These sizes refer to the annotation-specific input sets rather than the initial disease-intersection totals. Oxidative-stress response (GO:0006979) included 9/78 genes for ZBM (adjusted P =5.98 × 10^-4^), 106/835 for XKC (6.35 × 10^-50^), and 149/1,582 for XYT (7.06 × 10^-57^). Calcium homeostasis and voltage-gated calcium channel complex annotations were also represented. For the latter term (GO:0005891), GeneRatio was 15/78, 15/839, and 19/1,587, respectively. The annotation profiles are shown in Figure 2.

![](20261024_03_support/20261024_Figure2.png){height=8.4in}

**Figure 2:** Selected GO annotations. Panels show biological process (BP), molecular function (MF), and cellular component (CC). Terms are the union of the three lowest-adjusted-P terms per group within each ontology. Bar length shows GeneRatio. Colours denote ZBM (orange), XKC (teal), and XYT (purple). Missing bars denote terms absent from the exported significant-result table. The selected terms and their underlying statistics are provided in the accompanying source table.

The saved KEGG tables contained 83, 221, and 239 pathways with adjusted P <0.05 for ZBM, XKC, and XYT, using 74, 724, and 1,346 annotated input genes, respectively. Calcium signalling, MAPK signalling, and apoptosis were represented in all three datasets. The thyroid cancer pathway was present for XKC, with 16/724 genes (adjusted P =2.56 × 10^-8^), and XYT, with 20/1,346 genes (1.16 × 10^-7^), but was absent from the exported significant ZBM table. Lipid and atherosclerosis was the lowest-adjusted-P KEGG term for XKC and XYT, whereas nicotine addiction ranked first for ZBM (Figure 3). Pathway labels describe shared annotations and do not imply treatment of the named diseases.

![](20261024_03_support/20261024_Figure3.png)

**Figure 3:** Selected KEGG pathways. Selection uses the union of the five lowest-adjusted-P pathways per group, plus hsa05216 (thyroid cancer). Bars show GeneRatio; colours and missing-bar interpretation follow Figure 2. This descriptive display does not test differences between groups. Full enrichment results and figure source rows accompany the manuscript.

**Paired PTC transcriptome analysis.** The final analysis included 24 patients with 28 primary-tumour blocks and 27 normal-thyroid blocks, averaged to 48 patient–tissue observations. Of 20,848 measured genes, 666 met the differential-expression criteria: 347 were higher and 319 lower in tumour tissue (Figure 4). XYT included 111 differentially expressed genes among 1,518 measured candidates. The corresponding odds ratios were 2.67 against all measured genes (adjusted P =6.08 × 10^-16^) and 1.80 against the disease-restricted background (adjusted P =5.31 × 10^-7^). XKC included 69/805 differentially expressed genes, with odds ratios of 3.05 and 2.08 (adjusted P =2.44 × 10^-13^ and 5.31 × 10^-7^). ZBM included 4/77 differentially expressed genes; neither test was significant (adjusted P =0.278 and 0.482). POLB showed a small decrease (log2FC =−0.070; adjusted P =0.00271) that fell below the required fold-change threshold. CACNA1H and CYP3A4 showed small, non-significant changes (log2FC =−0.053 and −0.031; adjusted P =0.273 and 0.338). Thus, none of these three genes met both criteria.

![](20261025_04_review_support/20261025_Figure4.png)

**Figure 4:** Paired PTC transcriptome analysis in GSE60542 (24 patients). (A) Paired-limma volcano plot for 20,848 measured genes. Dashed lines mark adjusted P =0.05 and |log2FC| =1. Candidate sets overlap: XYT points (purple) were drawn first, followed by XKC (teal) and ZBM (orange), so the visible colour does not imply exclusive set membership. (B) Conditional odds ratios with two-sided exact 95% confidence intervals for over-representation among differentially expressed genes. Backgrounds comprise all measured genes or measured genes in the archived disease list. FDR labels show Benjamini–Hochberg-adjusted P values from the six one-sided Fisher tests.

β-Sitosterol (CID 222284) was represented in all three ingredient sets and was selected for further structural exploration. The same CID–receptor combinations occurred under different herb labels, reflecting the shared compound rather than chemically distinct herb-specific ligands. These repeated combinations were therefore not interpreted as evidence of differences in affinity between the herbal groups.

Supplementary smina outputs gave minimum scores of −9.12 kcal/mol for CYP3A4 (6UNI) and −9.68 kcal/mol for CaV3.2 (9AYG) (Figure 5). The 4P5A output was excluded from human-target interpretation. Receptor preparation removed the CYP3A4 haem cofactor, and the automatic box reference included both haem and the crystallographic inhibitor. For 9AYG, the reference contained multiple non-protein residue types without a verified antagonist-defined site. These limitations precluded interpreting the poses or scores as validated binding evidence. No reference-ligand redocking assessment or direct target-binding experiment was available.

![](20261025_04_review_support/20261025_Figure5.png)

**Figure 5:** Exploratory docking poses of β-sitosterol (CID 222284). (A) Human CYP3A4 (PDB 6UNI); (B) human CaV3.2/CACNA1H (PDB 9AYG). Left panels show the local view of the docked ligand (green sticks) and surrounding residues; right panels show the global receptor view. Poses correspond to the minimum supplementary smina scores of −9.12 and −9.68 kcal/mol, respectively. The bacterial 4P5A entry was excluded from human-target interpretation (see Methodology). These poses are exploratory structural outputs and are not validated binding evidence.

# Discussion

This study linked the selected XYT core herbs to thyroid-disease-associated genes and examined those genes in paired PTC tissue. The functional annotations involved oxidative stress, calcium regulation, and signalling. In GSE60542, XYT and XKC candidates were over-represented among genes altered in tumours, including when the comparison was restricted to measured genes already present in the disease list. The smaller odds ratios under that restricted background show that part of the enrichment depends on the choice of comparator genes. Persistence of the association supports a relationship to the PTC expression profile, while leaving the effects of herbal treatment untested.

The three gene sets are nested: both individual-herb intersections are contained in XYT. Broader coverage for XYT is therefore expected and cannot establish synergy. Although ZBM enrichment was not statistically significant, it included only 77 measured genes and four differentially expressed genes, giving a wide confidence interval. Its non-significant result does not establish lack of activity or a statistically demonstrated difference from XYT or XKC. Nor do the different overlap percentages or GeneRatios measure treatment effects.

The oxidative-stress annotations are consistent with the established relevance of redox regulation in thyroid cancer.^18^ Calcium-related annotations also warrant investigation, given experimental evidence linking TRPC1 to thyroid cancer cell migration and proliferation.^19^ However, TRPC1 and CACNA1H encode different channel types, and a shared calcium-related annotation does not connect them mechanistically. The current transcriptome analysis assessed whole candidate sets, not independent activity of each GO or KEGG pathway. Likewise, enrichment of the KEGG thyroid cancer pathway is partly conditioned on the use of a thyroid-tumour gene list to select the inputs.

The expression results did not provide strong gene-level support for the three preselected genes considered alongside docking. POLB was statistically different but changed by less than 0.1 log2 units; CACNA1H and CYP3A4 were not significant. Transcript abundance alone does not determine protein activity, but these findings do not justify promoting these genes as established therapeutic targets. CYP3A4 is involved in vandetanib metabolism,^20^ so a putative interaction may relate to drug metabolism rather than a tumour-selective effect. The archived docking results cannot resolve this question because CYP3A4 was prepared without its haem cofactor and the automatic box references did not establish validated binding sites. The same compound–receptor pair also remains the same chemical interaction regardless of which herb supplied the compound record.

Previous β-sitosterol studies provide useful context but concern different experimental systems. Cao and colleagues studied non-small-cell lung cancer,^21^ whereas Chen and colleagues examined Prunella vulgaris components in inflammatory thyrocyte models relevant to Hashimoto's disease.^22^ Neither study establishes the proposed interactions in PTC. Wnt and Keap1/Nrf2 signalling are relevant to thyroid cancer biology,^23,24^ but their modulation by the selected herbs was not measured here.

Several methodological limitations affect interpretation. The compound–gene table pools recorded and predicted associations, and its inclusion values are not comparable experimental confidence scores. BindingDB records were not restricted by organism or affinity, and a field-name inconsistency excluded TCMSP target records. These features may alter the retained target profile. The broad disease list has an unresolved historical search syntax. GO and KEGG analyses used annotation-wide backgrounds rather than the genes eligible for selection from the compound–target resource; the disease-restricted sensitivity analysis addressed transcriptome overlap, not this separate functional-enrichment limitation. Correlation among genes and shared database annotation can also affect gene-set tests. The external expression evidence comes from one microarray cohort and applies to PTC, not all thyroid cancer subtypes. Differences in tissue composition may contribute to the observed expression changes. There were no treatment-exposed samples, measurements of constituent concentrations, or direct binding and functional assays.

The next experimental step is to test chemically characterised preparations or authenticated compounds in thyroid cancer models with non-malignant controls, followed by target perturbation and binding or functional assays. Any renewed docking analysis should first use appropriate human structures, cofactors, defined binding sites, and reference-ligand controls.

# Conclusion

XYT and XKC candidate genes were over-represented among genes altered in paired PTC tissue, including against a disease-restricted background. ZBM enrichment was inconclusive, and the three preselected genes did not meet the combined differential-expression criteria. The findings narrow the hypotheses for further testing; they do not establish herb efficacy, synergy, or the proposed compound–protein interactions.

# Conflict of Interest

Zhiyu Luo is an employee of Jiangsu Sanshu Biotechnology Co., Ltd. The remaining authors declare no conflicts of interest related to this work.

# Funding

This work was supported by the construction project of the inheritance studio of Wang Wenchun, a nationally renowned expert in traditional Chinese medicine, under the National Administration of Traditional Chinese Medicine (Document No. Guo Zhong Yi Yao Ren Jiao Fa [2010] 59); and Shanghai Clinical Research Center of Traditional Chinese Medicine Oncology Science and Technology Commission of Shanghai Municipality (21MC1930500). Additional support was provided by the 2025 Research Project of the Shanghai Research Center for the Theory of Education for Older Adults (ZX0551-25) and the 2025 project “Policies and Pathways for Training Interdisciplinary Digital Health Professionals in Shanghai in the Context of the Digital and Intelligent Transformation of Healthcare” (2025JZ02). The listed funding sources provided support for this work.

# Patient Consent

No new consent was obtained for this secondary analysis of publicly available, de-identified data. No participants were recruited or contacted.

# Ethics Approval

This study used database records, molecular structures, and publicly available, de-identified GEO expression data. It involved no new participant recruitment, intervention, or animal experiments.

# Author Contributions

Tao Li: Data collection; data curation; formal analysis. Jiacai Li: Supervision; writing – review and editing. Yajie Ding: Data curation; writing – review and editing. Zhiyu Luo: Formal analysis. Xing Zhang: Data curation; writing – review and editing. Pengfei Yang: Supervision. All authors have reviewed and approved the final manuscript.

# Data Availability

Gene lists, sample and pairing manifests, differential-expression results, enrichment tests, and quality-control figures are included in the supplementary material. Expression data and platform annotations are available from GEO under GSE60542 and GPL570. The underlying network-pharmacology data are available from the corresponding author upon reasonable request. Analysis code is not publicly available.

# References

1. Cabanillas ME, McFadden DG, Durante C. Thyroid cancer. Lancet. 2016;388(10061):2783–95. doi:10.1016/S0140-6736(16)30172-6.

2. Haugen BR, Alexander EK, Bible KC, Doherty GM, Mandel SJ, Nikiforov YE, et al. 2015 American Thyroid Association Management Guidelines for Adult Patients with Thyroid Nodules and Differentiated Thyroid Cancer. Thyroid. 2016;26(1):1–133. doi:10.1089/thy.2015.0020.

3. Xing M. Molecular pathogenesis and mechanisms of thyroid cancer. Nat Rev Cancer. 2013;13:184–99. doi:10.1038/nrc3431.

4. Ru J, Li P, Wang J, Zhou W, Li B, Huang C, et al. TCMSP: a database of systems pharmacology for drug discovery from herbal medicines. J Cheminform. 2014;6:13. doi:10.1186/1758-2946-6-13.

5. Kim S, Chen J, Cheng T, Gindulyte A, He J, He S, et al. PubChem 2023 update. Nucleic Acids Res. 2023;51(D1):D1373–80. doi:10.1093/nar/gkac956.

6. Gilson MK, Liu T, Baitaluk M, Nicola G, Hwang L, Chong J. BindingDB in 2015: a public database for medicinal chemistry, computational chemistry and systems pharmacology. Nucleic Acids Res. 2016;44(D1):D1045–53. doi:10.1093/nar/gkv1072.

7. Kong X, Liu C, Zhang Z, Cheng M, Mei Z, Li X, et al. BATMAN-TCM 2.0: an enhanced integrative database for known and predicted interactions between traditional Chinese medicine ingredients and target proteins. Nucleic Acids Res. 2024;52(D1):D1110–20. doi:10.1093/nar/gkad926.

8. Harding SD, Armstrong JF, Faccenda E, Southan C, Alexander SPH, Davenport AP, et al. The IUPHAR/BPS Guide to PHARMACOLOGY in 2024. Nucleic Acids Res. 2024;52(D1):D1438–49. doi:10.1093/nar/gkad944.

9. Safran M, Dalah I, Alexander J, Rosen N, Iny Stein T, Shmoish M, et al. GeneCards Version 3: the human gene integrator. Database (Oxford). 2010;2010:baq020. doi:10.1093/database/baq020.

10. Tarabichi M, Saiselet M, Trésallet C, Hoang C, Larsimont D, Andry G, et al. Revisiting the transcriptional analysis of primary tumours and associated nodal metastases with enhanced biological and statistical controls: application to thyroid cancer. Br J Cancer. 2015;112(10):1665–74. doi:10.1038/bjc.2014.665.

11. Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK. limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Res. 2015;43(7):e47. doi:10.1093/nar/gkv007.

12. Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS. 2012;16(5):284–7. doi:10.1089/omi.2011.0118.

13. O'Boyle NM, Banck M, James CA, Morley C, Vandermeersch T, Hutchison GR. Open Babel: an open chemical toolbox. J Cheminform. 2011;3:33. doi:10.1186/1758-2946-3-33.

14. Burley SK, Bhikadiya C, Bi C, Bittrich S, Chen L, Crichlow GV, et al. RCSB Protein Data Bank: powerful new tools for exploring 3D structures of biological macromolecules for basic and applied research and education in fundamental biology, biomedicine, biotechnology, bioengineering and energy sciences. Nucleic Acids Res. 2021;49(D1):D437–51. doi:10.1093/nar/gkaa1038.

15. Trott O, Olson AJ. AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. J Comput Chem. 2010;31(2):455–61. doi:10.1002/jcc.21334.

16. Koes DR, Baumgartner MP, Camacho CJ. Lessons learned in empirical scoring with smina from the CSAR 2011 benchmarking exercise. J Chem Inf Model. 2013;53(8):1893–904. doi:10.1021/ci300604z.

17. Huang J, Fan X, Jin X, Lyu C, Guo Q, Liu T, et al. Structural basis for human CaV3.2 inhibition by selective antagonists. Cell Res. 2024;34:440–50. doi:10.1038/s41422-024-00959-8.

18. Cazarin J, Dupuy C, Pires de Carvalho D. Redox homeostasis in thyroid cancer: implications in Na+/I− symporter (NIS) regulation. Int J Mol Sci. 2022;23(11):6129. doi:10.3390/ijms23116129.

19. Asghar MY, Magnusson M, Kemppainen K, Sukumaran P, Löf C, Pulli I, et al. Transient receptor potential canonical 1 (TRPC1) channels as regulators of sphingolipid and VEGF receptor expression: implications for thyroid cancer cell migration and proliferation. J Biol Chem. 2015;290(26):16116–31. doi:10.1074/jbc.M115.643668.

20. Indra R, Pompach P, Martinek V, Takácsová P, Vavrová K, Heger Z, et al. Identification of human enzymes oxidizing the anti-thyroid-cancer drug vandetanib and explanation of the high efficiency of cytochrome P450 3A4 in its oxidation. Int J Mol Sci. 2019;20(14):3392. doi:10.3390/ijms20143392.

21. Cao W, Yuan F, Liu T, Yin R. Network pharmacology analysis, molecular docking integrated experimental verification reveal β-sitosterol as the active anti-NSCLC ingredient of Polygonatum cyrtonema Hua by suppression of PI3K/Akt/HIF-1α signaling pathway. J Ethnopharmacol. 2024;328:117900. doi:10.1016/j.jep.2024.117900.

22. Chen Y, Jiang B, Qu C, Jiang C, Zhang C, Wang Y, et al. Bioactive components in Prunella vulgaris for treating Hashimoto's disease via regulation of innate immune response in human thyrocytes. Heliyon. 2024;10(16):e36103. doi:10.1016/j.heliyon.2024.e36103.

23. Sastre-Perona A, Santisteban P. Role of the Wnt pathway in thyroid cancer. Front Endocrinol (Lausanne). 2012;3:31. doi:10.3389/fendo.2012.00031.

24. Renaud CO, Ziros PG, Chartoumpekis DV, Bongiovanni M, Sykiotis GP. Keap1/Nrf2 signaling: a new player in thyroid pathophysiology and thyroid cancer. Front Endocrinol (Lausanne). 2019;10:510. doi:10.3389/fendo.2019.00510.
