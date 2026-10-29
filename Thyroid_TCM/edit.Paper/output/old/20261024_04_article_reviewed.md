# Network pharmacology of core herbs selected from Xiao Ying Tang: thyroid-cancer-associated targets and functional enrichment

Tao Li^1^\*, Yajie Ding^2^\*, Zhiyu Luo^3^\*, Xing Zhang^4^, and Pengfei Yang^5^

^1^[AUTHOR CONFIRM: affiliation for Tao Li]

^2^[AUTHOR CONFIRM: affiliation for Yajie Ding]

^3^[AUTHOR CONFIRM: affiliation for Zhiyu Luo]

^4^[AUTHOR CONFIRM: affiliation for Xing Zhang]

^5^Clinical College of Chinese Medicine, Gansu University of Chinese Medicine, No. 35 Dingxi East Road, Chengguan District, Lanzhou City, Gansu Province, People's Republic of China

\*These authors contributed equally to this work.

Correspondence: Pengfei Yang, Clinical College of Chinese Medicine, Gansu University of Chinese Medicine, No. 35 Dingxi East Road, Chengguan District, Lanzhou City, Gansu Province, People's Republic of China. Tel: +8618293195696; Email: ypf@gszy.edu.cn

**REVIEWED WORKING DRAFT — 24 October 2026, revision 2.** The study investigates core herbs selected from Xiao Ying Tang for thyroid disease. XYT denotes this selected herb set throughout. Yellow-highlighted author queries concern reporting details still to be completed. The revised title emphasises the functional analysis; the unvalidated ancillary docking is reported with its limitations.

# Abstract

**Objective:** To explore thyroid-cancer-associated targets and pathways of core herbs selected from Xiao Ying Tang (XYT), with Zhe Bei Mu (ZBM) and Xia Ku Cao (XKC) as individual-herb comparators.

**Study design:** Descriptive in silico network pharmacology study with ancillary docking assessment.

**Place & duration of study:** [AUTHOR CONFIRM: institution and actual study dates, including the supplementary docking stage.]

**Methodology:** TCMSP-derived ingredients were linked to compound-associated genes and intersected with a GeneCards-derived disease set. Gene Ontology (GO) and Kyoto Encyclopedia of Genes and Genomes (KEGG) over-representation analyses used Benjamini–Hochberg adjustment. Functional profiles were compared descriptively; compound coverage, annotated backgrounds, and receptor preparation were checked.

**Results:** Among 12,987 disease-associated genes, 78, 842, and 1,595 overlapped with the 107, 1,094, and 2,131 ingredient-associated genes of ZBM, XKC, and XYT, respectively. Oxidative-stress response and calcium-related annotations were enriched. The KEGG thyroid cancer pathway included 16/724 annotated genes in XKC (adjusted P =2.56 × 10^-8^) and 20/1,346 in XYT (1.16 × 10^-7^). The ZBM and XKC disease-overlapping gene sets were contained within XYT. Docking outputs did not provide validated binding evidence.

**Conclusion:** The selected XYT core herbs showed broader gene coverage than the individual-herb comparators and shared annotations involving oxidative stress, calcium regulation, and signalling. The findings identify candidate processes for experimental study; the nested gene sets do not establish comparative efficacy.

**Keywords:** Thyroid Neoplasms; Medicine, Chinese Traditional; Gene Ontology; Computational Biology

# Introduction

Thyroid cancer comprises biologically distinct malignancies with different clinical courses and treatment requirements.^1^ Surgery, selective use of radioactive iodine, and appropriate thyroid hormone management are established components of differentiated thyroid cancer care.^2^ Persistent, metastatic, or treatment-resistant disease nevertheless motivates investigation of additional therapeutic mechanisms. Molecular alterations and signalling dependencies differ between thyroid cancer subtypes; therefore, a general thyroid-cancer gene list cannot establish a mechanism shared by all subtypes.^3^

Traditional Chinese medicine (TCM) ingredients can be investigated computationally to generate candidate compound–target hypotheses. However, a database association does not demonstrate that an ingredient is present at an active concentration in a preparation, reaches a tumour, or modifies a target in vivo. Enrichment and docking results consequently require careful separation from evidence of treatment efficacy.

Xiao Ying Tang provides the prescription context for this investigation. The research focused on core herbs selected from the prescription for thyroid disease, with the abbreviation XYT referring to this selected herb set throughout. Zhe Bei Mu (ZBM; Fritillariae Thunbergii Bulbus) and Xia Ku Cao (XKC; Prunellae Spica) were also analysed individually to describe their individual annotation profiles alongside the broader selection. The objective was to identify thyroid-cancer-associated gene overlap and functional themes across these ingredient sets and to examine candidate compound–protein interactions through exploratory docking.

# Methodology

An in silico study integrated herbal ingredient records, compound–gene associations, disease-associated genes, functional enrichment, and molecular docking. The selected core-herb set was the unit of analysis for XYT; ZBM and XKC were individual-herb comparators. No patients, animals, or clinical outcomes were studied. [AUTHOR CONFIRM: research institution, study period, database retrieval dates, and software/database versions.]

Ingredient records were obtained from TCMSP and matched to PubChem compound identifiers (CIDs).^4,5^ Lower-case ingredient names were joined to a compound-synonym table and then to an integrated compound–gene table. Records satisfying the source-table criterion `Gene_Value > 80` were retained, and gene symbols were aggregated by CID. BindingDB was among the interaction resources described for target collection.^6^ [AUTHOR CONFIRM: complete source list, the definition and scale of Gene_Value, and the experimental or predictive evidence underlying the associations.] No oral-bioavailability or drug-likeness cut-off was applied in the available analysis script.

The core herbs were selected by the investigators from the prescription context for their relevance to thyroid disease. The analysed XYT ingredient table contained the herb labels E Zhu, Bai Shao, Dang Gui, Chai Hu, Qing Pi, Huang Yao Zi, Xia Ku Cao, Zao Jiao Ci, Mu Dan Pi, Bai Zhu, Fu Ling, Bai Mao Xia Ku Cao, Xuan Shen, and Zhe Bei Mu. This inventory defines the computational exposure evaluated here. ZBM and XKC were extracted as individual-herb datasets using their Chinese herb names. All three sets were processed using the same compound-matching and gene-intersection workflow. Distinct CIDs were counted at three stages: all mapped ingredients, ingredients with at least one retained gene association, and ingredients associated with at least one disease-overlapping gene. Duplicate herb–compound records were collapsed within each set for these counts. The resulting database coverage did not encode herb doses or relative constituent abundance.

Disease-associated genes were obtained from a GeneCards export containing 12,987 unique gene symbols.^7^ The same disease set was used for all three comparisons, with no additional relevance-score filter in the intersection step. [AUTHOR CONFIRM: exact GeneCards query and any filtering applied before export.] Ingredient-associated symbols were split on semicolons and deduplicated before intersection. Set counts and membership were verified against the saved common-gene lists. Disease-only, ingredient-only, and shared genes were recorded separately; an association with thyroid disease was not taken as proof of causal involvement or tumour-specific expression.

Gene symbols were mapped to Entrez identifiers with org.Hs.eg.db, followed by GO and KEGG over-representation analysis using clusterProfiler in R.^8^ GO analysis covered biological process (BP), molecular function (MF), and cellular component (CC) using `ont = "ALL"`; KEGG analysis used the human annotation (`organism = "hsa"`). Both analyses specified Benjamini–Hochberg adjustment, `pvalueCutoff = 0.05`, and `qvalueCutoff = 0.2`. No custom background universe was supplied. The exported BgRatio fields indicated effective backgrounds of 18,888 genes for BP, 18,522 for MF, 19,894 for CC, and 8,859 for KEGG. Over-representation used an upper-tail hypergeometric test. Reported terms had adjusted P <0.05. As a numerical check, raw P values were independently recalculated from the exported Count, GeneRatio, and BgRatio fields. The saved adjusted P values were retained because the complete pre-filtering hypothesis lists were unavailable for independent reconstruction of the BH adjustment. [AUTHOR CONFIRM: annotation versions and the handling of unmapped or one-to-many identifiers.]

For descriptive visualisation, GO terms were selected as the union of the three lowest-adjusted-P terms for each dataset within each ontology. KEGG pathways used the five lowest-adjusted-P pathways per dataset, with the thyroid cancer pathway (hsa05216) additionally included because of the study question. GeneRatio was calculated from the numerator and denominator recorded in each result row. The denominator represents the annotated input genes used for that analysis and can differ between ontologies. Missing bars indicate absence from the exported significant-result table, not absence of biological involvement. No between-dataset test of enrichment strength or therapeutic efficacy was performed. Nominally overlapping terms were not treated as independent biological mechanisms, and a term significant in one set but absent from another significant-result table was not interpreted as a statistically significant between-set difference.

Ligand structures were converted using Open Babel, and receptor structures were retrieved from RCSB PDB.^9,10^ The initial screening used AutoDock Vina; supplementary docking of β-sitosterol used smina, a Vina-derived program.^11,12^ Receptor identity and preparation were checked before interpretation. PDB 6UNI corresponds to human CYP3A4, and 9AYG to human CaV3.2/CACNA1H in an apo state.^10,13^ The 4P5A entry was excluded because it represents bacterial ThyX/PolB rather than human DNA polymerase beta. The supplementary protocol lacked a complete recoverable command for the search box, random seed, exhaustiveness, and scoring settings. Consequently, scores were treated as exploratory outputs; no universal score threshold was used to classify binding strength or establish inhibition.

# Results

The ingredient tables contained 16, 68, and 803 distinct non-missing PubChem CIDs for ZBM, XKC, and XYT, respectively. Of these, 8, 57, and 501 had retained gene associations; 8, 57, and 492 were associated with at least one gene in the disease intersection. These counts distinguish database ingredient coverage from the compounds that contributed disease-associated genes. The corresponding unique ingredient-associated gene totals were 107, 1,094, and 2,131. Of these, 78/107 (72.9%), 842/1,094 (77.0%), and 1,595/2,131 (74.8%) overlapped with the common disease set. Disease-only counts were 12,909, 12,145, and 11,392, while ingredient-only counts were 29, 252, and 536 (Figure 1). All 78 ZBM and all 842 XKC disease-overlapping genes were contained in the 1,595-gene XYT intersection. XYT therefore had the largest absolute coverage in a nested comparison; the overlap percentages describe these sets without testing differences between them.

![](20261024_03_support/20261024_Figure1.png)

**Figure 1:** Gene-set intersections for ZBM, XKC, and the XYT core-herb set. Shared disease genes total 78, 842, and 1,595, respectively. All comparisons use the same 12,987-gene disease set. Circle areas are schematic and are not proportional to set sizes. XYT denotes the core herbs selected from Xiao Ying Tang for thyroid disease.

The saved GO result tables contained 1,026, 3,925, and 4,782 terms with adjusted P <0.05 for ZBM, XKC, and XYT. Their annotated input sizes were 78/78/78 for ZBM, 835/840/839 for XKC, and 1,582/1,588/1,587 for XYT in BP/MF/CC order. These sizes refer to the annotation-specific input sets rather than the initial disease-intersection totals. Oxidative-stress response (GO:0006979) included 9/78 genes for ZBM (adjusted P =5.98 × 10^-4^), 106/835 for XKC (6.35 × 10^-50^), and 149/1,582 for XYT (7.06 × 10^-57^). Calcium homeostasis and voltage-gated calcium channel complex annotations were also represented. For the latter term (GO:0005891), GeneRatio was 15/78, 15/839, and 19/1,587, respectively. These results illustrate why absolute gene counts and within-set proportions must be distinguished (Figure 2).

![](20261024_03_support/20261024_Figure2.png){height=8.4in}

**Figure 2:** Selected GO annotations. Panels show biological process (BP), molecular function (MF), and cellular component (CC). Terms are the union of the three lowest-adjusted-P terms per group within each ontology. Bar length shows GeneRatio. Colours denote ZBM (orange), XKC (teal), and XYT (purple). Missing bars denote terms absent from the exported significant-result table. The selected terms and their underlying statistics are provided in the accompanying source table.

The saved KEGG tables contained 83, 221, and 239 pathways with adjusted P <0.05 for ZBM, XKC, and XYT, using 74, 724, and 1,346 annotated input genes, respectively. Calcium signalling, MAPK signalling, and apoptosis were represented in all three datasets. The thyroid cancer pathway was present for XKC, with 16/724 genes (adjusted P =2.56 × 10^-8^), and XYT, with 20/1,346 genes (1.16 × 10^-7^), but was absent from the exported significant ZBM table. Lipid and atherosclerosis was the lowest-adjusted-P KEGG term for XKC and XYT, whereas nicotine addiction ranked first for ZBM (Figure 3). Pathway labels describe shared annotations and do not imply treatment of the named diseases. Across the 10,276 exported GO and KEGG rows, raw hypergeometric P values were reproduced with a maximum absolute difference below 2 × 10^-16^; this check confirms the numerical consistency of the saved rows rather than independently validating the input associations.

![](20261024_03_support/20261024_Figure3.png)

**Figure 3:** Selected KEGG pathways. Selection uses the union of the five lowest-adjusted-P pathways per group, plus hsa05216 (thyroid cancer). Bars show GeneRatio; colours and missing-bar interpretation follow Figure 2. This descriptive display does not test differences between groups. Full enrichment results and figure source rows accompany the manuscript.

β-Sitosterol (CID 222284) was represented in all three ingredient sets and was selected for further structural exploration. The same CID–receptor combinations occurred under different herb labels, reflecting the shared compound rather than chemically distinct herb-specific ligands. These repeated combinations were therefore not interpreted as evidence of differences in affinity between the herbal groups.

Supplementary smina outputs gave minimum scores of −9.12 kcal/mol for CYP3A4 (6UNI) and −9.68 kcal/mol for CaV3.2 (9AYG). The 4P5A output was excluded from human-target interpretation. Receptor preparation removed the CYP3A4 haem cofactor, and the automatic box reference included both haem and the crystallographic inhibitor. For 9AYG, the reference contained multiple non-protein residue types without a verified antagonist-defined site. These limitations precluded interpreting the poses or scores as validated binding evidence. No reference-ligand redocking assessment or direct target-binding experiment was available.

# Discussion

The selected XYT core herbs covered a broad set of thyroid-cancer-associated genes, with enrichment involving oxidative-stress responses, calcium regulation, and signalling. The individual-herb profiles overlapped with these themes, indicating common annotated processes within the broader core-herb selection. The contribution of this comparison is to identify processes and candidate associations that can be taken forward for targeted testing; enrichment alone does not determine whether a pathway is activated or inhibited.

The distinction between absolute coverage and relative representation is central to interpreting the comparison. XYT had the largest number of disease-overlapping genes, but the overlap proportions were 74.8% for XYT, 77.0% for XKC, and 72.9% for ZBM. Similarly, the voltage-gated calcium channel complex annotation involved more genes in XYT than ZBM but a smaller GeneRatio. Both measures are informative: counts describe coverage, whereas GeneRatio describes the proportion of annotated input genes assigned to a term. Neither is a measure of treatment effect. The gene-level check showed that both individual-herb disease intersections were subsets of the XYT intersection. Their shared biological themes are therefore partly expected from the selection structure. The comparison describes coverage and annotation profiles, but does not isolate each herb’s marginal contribution, demonstrate independent treatment effects, or establish pharmacological synergy.

Redox regulation provides a biologically relevant context for the oxidative-stress annotations.^14^ Calcium regulation is also of interest because experimental work has implicated TRPC1 in thyroid cancer cell migration and proliferation.^15^ However, the GO enrichment does not identify the direction of regulation, and TRPC channels are distinct from the CaV3.2 channel encoded by CACNA1H. Linking these findings requires measurement of the relevant proteins, calcium responses, and downstream phenotypes in an appropriate thyroid cancer model. The KEGG thyroid cancer annotation in XKC and XYT identifies gene subsets for further investigation. Because the inputs were already selected by overlap with a thyroid-tumour gene list, enrichment of a thyroid-cancer-labelled pathway is not an independent validation of disease specificity. Broad or non-thyroid disease annotations can reflect shared signalling components.

CYP3A4 also requires a different biological interpretation from a tumour-selective therapeutic target. Its role in vandetanib metabolism has been experimentally examined.^16^ A hypothetical interaction with CYP3A4 could concern metabolism or drug–herb interactions; it does not by itself indicate anticancer benefit. Deleting its haem cofactor further weakens interpretation of the archived docking pose. Likewise, a negative docking score is an empirical scoring output, not proof of spontaneous binding, target inhibition, or clinical benefit. Differences between repeated runs of an identical ligand–receptor pair cannot be attributed to the herb from which that ligand was catalogued.

Published experiments on β-sitosterol provide context but do not validate this study's proposed targets. Cao and colleagues investigated β-sitosterol in non-small-cell lung cancer models, not thyroid cancer.^17^ Chen and colleagues studied Prunella vulgaris components in inflammatory thyrocyte models relevant to Hashimoto's disease; these are not thyroid cancer models, and β-sitosterol showed no detectable affinity for HSP90AA1 in their reported binding assay.^18^ That result does not test CYP3A4 or CACNA1H, but it illustrates why network predictions require direct confirmation. Wnt and Keap1/Nrf2 signalling provide broader thyroid-cancer context,^19,20^ without establishing that these pathways were modulated in the present study.

Experimental follow-up should use chemically characterised preparations or authenticated compounds, concentration–response testing in thyroid cancer models with non-malignant controls, and target perturbation alongside binding or functional assays. Docking should first use verified human receptor identities, appropriate cofactors, explicit binding sites, and reference compounds.

Several limitations should be considered. Compound–gene associations depend on database coverage, source evidence, and the confidence criterion used to retain records. The disease set is broad and lacks subtype-specific or expression-based validation. Over-representation results depend on annotation backgrounds, identifier mapping, and correlated terms. The annotation-wide backgrounds do not condition on prior selection of thyroid-tumour-associated genes or on which genes were accessible to the compound–target resource. Analyses using a justified eligible-gene background and sensitivity to the disease-gene selection would help distinguish ingredient-associated signals from database and selection effects. The present analysis did not perform these additional tests, and the nested ingredient-derived gene sets limit comparative interpretation. The study did not measure constituent concentrations, bioavailability, target engagement, cellular responses, or clinical outcomes. Docking additionally requires verified human receptor identity, appropriate cofactor treatment, explicit binding-site definition, and controlled validation. These limitations define the next steps without negating the descriptive network-pharmacology results.

# Conclusion

The selected XYT core herbs yielded a larger disease-associated gene set that contained the individual-herb intersections. Functional annotation highlighted oxidative stress, calcium regulation, and signalling as candidate themes. β-Sitosterol was a common candidate ingredient, but direct binding and antitumour activity were not established. The findings provide a basis for prioritising experimental studies of the selected core herbs and their candidate targets.

# Conflict of Interest

None to declare.

# Funding

This work was supported by the construction project of the inheritance studio of Wang Wenchun, a nationally renowned expert in traditional Chinese medicine, under the National Administration of Traditional Chinese Medicine (Document No. Guo Zhong Yi Yao Ren Jiao Fa [2010] 59); and Shanghai Clinical Research Center of Traditional Chinese Medicine Oncology Science and Technology Commission of Shanghai Municipality (21MC1930500). [AUTHOR CONFIRM: official award wording and funder roles.]

# Patient Consent

Not applicable: the reported analysis did not recruit patients or use individual patient data.

# Ethics Approval

The reported work analysed database records and molecular structures without recruiting human participants or using animals. [AUTHOR CONFIRM: institutional determination regarding ethics review or exemption and any documentation required by the journal.]

# Disclosure

[AUTHOR CONFIRM: relevant non-financial disclosures and accurate disclosure of computational/editorial assistance according to the intended journal's policy.]

# Author Contributions

[AUTHOR CONFIRM: provide individual contributions for Tao Li, Yajie Ding, Zhiyu Luo, Xing Zhang, and Pengfei Yang. Each author must review and approve the revised manuscript and accept accountability.]

# Data Availability

The accompanying analysis package contains gene-set counts, the herb inventory, enrichment tables, figure source data, and a reproducible checking script. Additional docking preparation checks are documented separately. [AUTHOR CONFIRM: public repository/accession or a justified data-access statement for submission, including the original compound–target source table and complete docking configurations.]

# References

1. Cabanillas ME, McFadden DG, Durante C. Thyroid cancer. Lancet. 2016;388(10061):2783–95. doi:10.1016/S0140-6736(16)30172-6.

2. Haugen BR, Alexander EK, Bible KC, Doherty GM, Mandel SJ, Nikiforov YE, et al. 2015 American Thyroid Association Management Guidelines for Adult Patients with Thyroid Nodules and Differentiated Thyroid Cancer. Thyroid. 2016;26(1):1–133. doi:10.1089/thy.2015.0020.

3. Xing M. Molecular pathogenesis and mechanisms of thyroid cancer. Nat Rev Cancer. 2013;13:184–99. doi:10.1038/nrc3431.

4. Ru J, Li P, Wang J, Zhou W, Li B, Huang C, et al. TCMSP: a database of systems pharmacology for drug discovery from herbal medicines. J Cheminform. 2014;6:13. doi:10.1186/1758-2946-6-13.

5. Kim S, Chen J, Cheng T, Gindulyte A, He J, He S, et al. PubChem 2023 update. Nucleic Acids Res. 2023;51(D1):D1373–80. doi:10.1093/nar/gkac956.

6. Gilson MK, Liu T, Baitaluk M, Nicola G, Hwang L, Chong J. BindingDB in 2015: a public database for medicinal chemistry, computational chemistry and systems pharmacology. Nucleic Acids Res. 2016;44(D1):D1045–53. doi:10.1093/nar/gkv1072.

7. Safran M, Dalah I, Alexander J, Rosen N, Iny Stein T, Shmoish M, et al. GeneCards Version 3: the human gene integrator. Database (Oxford). 2010;2010:baq020. doi:10.1093/database/baq020.

8. Yu G, Wang LG, Han Y, He QY. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS. 2012;16(5):284–7. doi:10.1089/omi.2011.0118.

9. O'Boyle NM, Banck M, James CA, Morley C, Vandermeersch T, Hutchison GR. Open Babel: an open chemical toolbox. J Cheminform. 2011;3:33. doi:10.1186/1758-2946-3-33.

10. Burley SK, Bhikadiya C, Bi C, Bittrich S, Chen L, Crichlow GV, et al. RCSB Protein Data Bank: powerful new tools for exploring 3D structures of biological macromolecules for basic and applied research and education in fundamental biology, biomedicine, biotechnology, bioengineering and energy sciences. Nucleic Acids Res. 2021;49(D1):D437–51. doi:10.1093/nar/gkaa1038.

11. Trott O, Olson AJ. AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. J Comput Chem. 2010;31(2):455–61. doi:10.1002/jcc.21334.

12. Koes DR, Baumgartner MP, Camacho CJ. Lessons learned in empirical scoring with smina from the CSAR 2011 benchmarking exercise. J Chem Inf Model. 2013;53(8):1893–904. doi:10.1021/ci300604z.

13. Huang J, Fan X, Jin X, Lyu C, Guo Q, Liu T, et al. Structural basis for human CaV3.2 inhibition by selective antagonists. Cell Res. 2024;34:440–50. doi:10.1038/s41422-024-00959-8.

14. Cazarin J, Dupuy C, Pires de Carvalho D. Redox homeostasis in thyroid cancer: implications in Na+/I− symporter (NIS) regulation. Int J Mol Sci. 2022;23(11):6129. doi:10.3390/ijms23116129.

15. Asghar MY, Magnusson M, Kemppainen K, Sukumaran P, Löf C, Pulli I, et al. Transient receptor potential canonical 1 (TRPC1) channels as regulators of sphingolipid and VEGF receptor expression: implications for thyroid cancer cell migration and proliferation. J Biol Chem. 2015;290(26):16116–31. doi:10.1074/jbc.M115.643668.

16. Indra R, Pompach P, Martinek V, Takácsová P, Vavrová K, Heger Z, et al. Identification of human enzymes oxidizing the anti-thyroid-cancer drug vandetanib and explanation of the high efficiency of cytochrome P450 3A4 in its oxidation. Int J Mol Sci. 2019;20(14):3392. doi:10.3390/ijms20143392.

17. Cao W, Yuan F, Liu T, Yin R. Network pharmacology analysis, molecular docking integrated experimental verification reveal β-sitosterol as the active anti-NSCLC ingredient of Polygonatum cyrtonema Hua by suppression of PI3K/Akt/HIF-1α signaling pathway. J Ethnopharmacol. 2024;328:117900. doi:10.1016/j.jep.2024.117900.

18. Chen Y, Jiang B, Qu C, Jiang C, Zhang C, Wang Y, et al. Bioactive components in Prunella vulgaris for treating Hashimoto's disease via regulation of innate immune response in human thyrocytes. Heliyon. 2024;10(16):e36103. doi:10.1016/j.heliyon.2024.e36103.

19. Sastre-Perona A, Santisteban P. Role of the Wnt pathway in thyroid cancer. Front Endocrinol (Lausanne). 2012;3:31. doi:10.3389/fendo.2012.00031.

20. Renaud CO, Ziros PG, Chartoumpekis DV, Bongiovanni M, Sykiotis GP. Keap1/Nrf2 signaling: a new player in thyroid pathophysiology and thyroid cancer. Front Endocrinol (Lausanne). 2019;10:510. doi:10.3389/fendo.2019.00510.
