# Exploratory target annotation of traditional Chinese medicine ingredient datasets in thyroid cancer: network pharmacology and a molecular docking assessment

Tao Li^1^\*, Yajie Ding^2^\*, Zhiyu Luo^3^\*, Xing Zhang^4^, and Pengfei Yang^5^

^1^[AUTHOR CONFIRM: affiliation for Tao Li]

^2^[AUTHOR CONFIRM: affiliation for Yajie Ding]

^3^[AUTHOR CONFIRM: affiliation for Zhiyu Luo]

^4^[AUTHOR CONFIRM: affiliation for Xing Zhang]

^5^Clinical College of Chinese Medicine, Gansu University of Chinese Medicine, No. 35 Dingxi East Road, Chengguan District, Lanzhou City, Gansu Province, People's Republic of China

\*These authors contributed equally to this work.

Correspondence: Pengfei Yang, Clinical College of Chinese Medicine, Gansu University of Chinese Medicine, No. 35 Dingxi East Road, Chengguan District, Lanzhou City, Gansu Province, People's Republic of China. Tel: +8618293195696; Email: ypf@gszy.edu.cn

**REVIEW WORKING DRAFT — 23 October 2026.** The supplied results have been checked against local scripts and tables. XYT remains a provisional dataset label because the analysed herb list differs from the stated Xiao Ying Tang prescription. The docking results do not establish target engagement. Author queries must be resolved before submission; this is not a submission-ready manuscript.

# Abstract

**Objective:** To characterise thyroid-cancer-associated gene annotations of Zhe Bei Mu (ZBM), Xia Ku Cao (XKC), and an archived multi-herb dataset labelled XYT, and assess the interpretability of their molecular docking results.

**Study design:** Exploratory computational study.

**Place & duration of study:** [AUTHOR CONFIRM: institution and actual study dates, including the supplementary docking stage.]

**Methodology:** Ingredient-associated genes were intersected with a GeneCards-derived disease set. Archived Gene Ontology (GO) and Kyoto Encyclopedia of Genes and Genomes (KEGG) results were reviewed using Benjamini–Hochberg-adjusted P values. Compound identities, herb selection, receptor annotations, and docking preparation records were checked.

**Results:** The common disease set contained 12,987 genes. ZBM, XKC, and XYT contained 107, 1,094, and 2,131 ingredient-associated genes, of which 78, 842, and 1,595 overlapped with the disease set. Enriched annotations included oxidative-stress responses, calcium regulation, and signalling pathways. Larger gene counts did not establish greater therapeutic activity. The analysed XYT herb inventory differed from the stated prescription. The proposed POLB receptor, PDB 4P5A, was a bacterial protein and was excluded. Supplementary docking used smina; receptor preparation removed the CYP3A4 haem cofactor, and binding-site definition remained insufficiently validated.

**Conclusion:** The datasets generated hypotheses about annotated biological processes but did not demonstrate anti-thyroid-cancer efficacy or superiority of XYT. Prescription reconciliation, target-source verification, and validated docking are required before compound–target mechanisms can be inferred.

**Keywords:** Thyroid Neoplasms; Medicine, Chinese Traditional; Molecular Docking Simulation; Gene Ontology; Computational Biology

# Introduction

Thyroid cancer comprises biologically distinct malignancies with different clinical courses and treatment requirements.^1^ Surgery, selective use of radioactive iodine, and appropriate thyroid hormone management are established components of differentiated thyroid cancer care.^2^ Persistent, metastatic, or treatment-resistant disease nevertheless motivates investigation of additional therapeutic mechanisms. Molecular alterations and signalling dependencies differ between thyroid cancer subtypes; therefore, a general thyroid-cancer gene list cannot establish a mechanism shared by all subtypes.^3^

Traditional Chinese medicine (TCM) ingredients can be investigated computationally to generate candidate compound–target hypotheses. However, a database association does not demonstrate that an ingredient is present at an active concentration in a preparation, reaches a tumour, or modifies a target in vivo. Enrichment and docking results consequently require careful separation from evidence of treatment efficacy.

The present study examined Zhe Bei Mu (ZBM; Fritillariae Thunbergii Bulbus), Xia Ku Cao (XKC; Prunellae Spica), and an archived multi-herb dataset labelled XYT. ZBM and XKC are individual herbal medicines, whereas XYT was intended to represent a prescription. Because the archived ingredient inventory and the manuscript prescription were inconsistent, XYT is retained here solely as a dataset identifier. The objective was to describe their overlap with thyroid-cancer-associated genes, assess functional annotations, and identify the limits of the available docking evidence.

# Methodology

This was an exploratory computational analysis of archived ingredient, gene, enrichment, and docking records. No patients, animals, or clinical outcomes were studied. [AUTHOR CONFIRM: original research institution, study period, database retrieval dates, and software/database versions.] Directory names were not treated as sufficient evidence of research dates.

Ingredient records originated from TCMSP, and compound identifiers were matched to PubChem CIDs.^4,5^ The archived script joined lower-case ingredient names to a compound-synonym table and then to a local compound–gene table, retaining records with `Gene_Value > 80`. The original methods named BindingDB and Guide to Pharmacology as sources; however, source-level provenance and the definition of `Gene_Value` were not recoverable from the inspected downstream tables. BindingDB is a relevant interaction resource,^6^ but the provenance of each retained association must be documented. The unverified database designation “MetaBat” was removed. The inspected script did not impose oral-bioavailability or drug-likeness cut-offs, and such filters were not assumed.

ZBM and XKC were selected by exact Chinese herb names. The multi-herb selection used a regular expression. Its output contained 14 herb labels: E Zhu, Bai Shao, Dang Gui, Chai Hu, Qing Pi, Huang Yao Zi, Xia Ku Cao, Zao Jiao Ci, Mu Dan Pi, Bai Zhu, Fu Ling, Bai Mao Xia Ku Cao, Xuan Shen, and Zhe Bei Mu. This list differs from the 11-herb Xiao Ying Tang prescription originally described. In particular, Bai Mao Xia Ku Cao was included alongside Xia Ku Cao. [AUTHOR CONFIRM: intended prescription, authenticated botanical identities, plant parts, preparation, and doses; repeat affected analyses if the selected inventory is incorrect.] No claim about the named clinical prescription was made from this unreconciled dataset.

The archived GeneCards export contained 12,987 unique gene symbols.^7^ Its filename indicated a thyroid-tumour search, but the exact query, retrieval date, and upstream relevance-score selection require confirmation. The intersection script used the exported set without a further relevance-score filter. Ingredient-associated symbols were split on semicolons and deduplicated. Set intersections were independently recalculated and checked against the saved common-gene files. Counts of disease-only genes, ingredient-only genes, and shared genes were distinguished explicitly. A disease association was not interpreted as evidence of tumour-specific expression, causal involvement, or druggability.

The archived R analysis used clusterProfiler and org.Hs.eg.db to map symbols to Entrez identifiers and perform GO and KEGG over-representation analyses.^8^ The script specified `ont = "ALL"` for GO, `organism = "hsa"` for KEGG, Benjamini–Hochberg adjustment, `pvalueCutoff = 0.05`, and `qvalueCutoff = 0.2`. No custom background universe was supplied. Interpretation was restricted to exported results with adjusted P <0.05. Enrichment tests were not rerun during revision; the values reported here are from the existing result files. Version-dependent annotation backgrounds, unmapped identifiers, and one-to-many mappings remain to be documented.

For descriptive visualisation, GO terms were selected as the union of the three lowest-adjusted-P terms for each dataset within each ontology. KEGG pathways used the five lowest-adjusted-P pathways per dataset, with the thyroid cancer pathway (hsa05216) additionally included because of the study question. GeneRatio was calculated from the numerator and denominator recorded in each result row. The denominator represents the annotated input genes used for that analysis and can differ between ontologies. Missing bars indicate absence from the exported significant-result table, not absence of biological involvement. No between-dataset test of enrichment strength or therapeutic efficacy was performed.

The original docking workflow used Open Babel for structure conversion, RCSB PDB structures, and AutoDock Vina.^9–11^ The supplementary archives instead identified smina, a Vina-derived docking program.^12^ Protein identity, organism, reference ligands, preparation logs, and saved scores were inspected. PDB identifiers were treated as structure accessions, not gene names. PDB 4P5A was excluded because it represents bacterial ThyX/PolB rather than human DNA polymerase beta. PDB 6UNI represents human CYP3A4, and 9AYG represents human CaV3.2/CACNA1H in an apo state.^10,13^ The supplementary archives lacked a complete recoverable command documenting search-box coordinates, random seed, exhaustiveness, and scoring settings. The original assertion of exhaustiveness 8 could therefore not be verified. These runs were assessed as unvalidated computational outputs rather than experimental confirmation.

# Results

The archived ingredient tables contained 16, 68, and 803 distinct non-missing PubChem CIDs for ZBM, XKC, and XYT, respectively. These are database identifiers, not experimentally measured constituents. The corresponding unique ingredient-associated gene totals were 107, 1,094, and 2,131. Of these, 78/107 (72.9%), 842/1,094 (77.0%), and 1,595/2,131 (74.8%) overlapped with the same disease set. The disease-only counts were 12,909, 12,145, and 11,392; the ingredient-only counts were 29, 252, and 536 (Figure 1). Thus, the values previously described as total genes were disjoint Venn-region counts.

![](20261023_03_support/20261023_Figure1.png)

**Figure 1:** Recalculated set intersections for ZBM, XKC, and the provisionally labelled XYT dataset. Shared disease genes total 78, 842, and 1,595, respectively. All comparisons use the same 12,987-gene disease set. Circle areas are schematic and are not proportional to set sizes. The asterisk identifies the unreconciled XYT herb inventory.

The exported GO tables contained 1,026, 3,925, and 4,782 significant terms for ZBM, XKC, and XYT. Oxidative-stress response (GO:0006979) included 9/78 genes for ZBM (adjusted P =5.98 × 10^-4^), 106/835 for XKC (6.35 × 10^-50^), and 149/1,582 for XYT (7.06 × 10^-57^). Calcium homeostasis and voltage-gated calcium channel complex annotations were also represented. For the latter term (GO:0005891), GeneRatio was 15/78, 15/839, and 19/1,587, respectively. These results illustrate why absolute gene counts and within-set proportions must be distinguished (Figure 2).

![](20261023_03_support/20261023_Figure2.png){height=8.4in}

**Figure 2:** Selected GO annotations from the archived significant-result tables. Panels show biological process (BP), molecular function (MF), and cellular component (CC). Terms are the union of the three lowest-adjusted-P terms per dataset within each ontology. Bar length shows GeneRatio, not effect size or therapeutic activity. Colours denote ZBM (orange), XKC (teal), and XYT (purple; provisional label). Missing bars denote terms absent from the exported significant-result table. Source rows are supplied in the accompanying audit files.

The exported KEGG tables contained 83, 221, and 239 significant pathways for ZBM, XKC, and XYT. Calcium signalling, MAPK signalling, and apoptosis were represented in all three datasets. The thyroid cancer pathway was present for XKC, with 16/724 genes (adjusted P =2.56 × 10^-8^), and XYT, with 20/1,346 genes (1.16 × 10^-7^), but was absent from the exported significant ZBM table. Lipid and atherosclerosis was the lowest-adjusted-P KEGG term for XKC and XYT, whereas nicotine addiction ranked first for ZBM (Figure 3). Pathway labels describe shared gene annotations and do not imply that the preparations treat the named diseases.

![](20261023_03_support/20261023_Figure3.png)

**Figure 3:** KEGG pathways selected from the archived results. Selection uses the union of the five lowest-adjusted-P pathways per dataset, plus hsa05216 (thyroid cancer). Bars show the recorded GeneRatio; colours and missing-bar interpretation follow Figure 2. This descriptive display does not test differences between datasets. Full exported results and figure source rows accompany the revised manuscript.

β-Sitosterol (CID 222284) was present in all three ingredient datasets. The original docking comparison therefore included repeated docking of the same chemical entity to the same structures under different herb labels, rather than independent formulation-specific ligands. The original “top 20 compounds” figure showed 17 labelled compound–structure pairs, all for CID 222284, and included the incorrect 4P5A–POLB assignment. It was removed from mechanistic interpretation.

The supplementary smina score files contained best scores of −9.23, −9.12, and −9.68 kcal/mol for 4P5A, 6UNI, and 9AYG. The 4P5A result was excluded on protein-identity grounds. For 6UNI, the log documented removal of heteroatoms, including haem; its autobox reference included both haem and the crystallographic inhibitor. For 9AYG, the reference file contained multiple non-protein residue types rather than a verified antagonist-defined binding site. The remaining scores were therefore retained only in the audit record, without claims of affinity, inhibition, or rank ordering between targets. No validated binding assay, reference-ligand redocking assessment, or biological efficacy experiment was available in the inspected study materials.

# Discussion

The principal finding was an overlap between database-derived ingredient associations and a broad thyroid-tumour gene set, with annotations involving calcium regulation, stress responses, and signalling. Such results can guide hypothesis generation, but they cannot establish pathway activation or inhibition. Redox regulation is relevant to thyroid physiology and cancer,^14^ and experimental work has implicated TRPC1 in thyroid cancer cell behaviour.^15^ Neither observation proves that the studied preparations inhibit calcium channels or suppress angiogenesis. In particular, TRPC channels and the CaV3.2 channel encoded by CACNA1H are distinct channel families.

The larger XYT gene list should not be interpreted as evidence of greater efficacy. It combines many more database compounds and includes both ZBM and XKC in its actual archived inventory. The comparisons are therefore overlapping, unequal ingredient sets rather than independent treatment groups. Shared ingredients and highly annotated compounds can increase both gene coverage and the number of significant terms. The unresolved prescription mismatch additionally prevents transfer of these results to the named Xiao Ying Tang preparation.

The receptor audit identified a substantive limitation of gene-name-only PDB retrieval. The local search used free-text gene queries without an explicit human-organism or UniProt-identity restriction, allowing bacterial PolB to be assigned to human POLB. Such mapping errors invalidate the corresponding mechanistic inference regardless of a favourable docking score. Correct mapping must precede any new docking or structure-based target prioritisation.

CYP3A4 also requires a different biological interpretation from a tumour-selective therapeutic target. Its role in vandetanib metabolism has been experimentally examined.^16^ A hypothetical interaction with CYP3A4 could concern metabolism or drug–herb interactions; it does not by itself indicate anticancer benefit. Deleting its haem cofactor further weakens interpretation of the archived docking pose. Likewise, a negative docking score is an empirical scoring output, not proof of spontaneous binding, target inhibition, or clinical benefit. Differences between repeated runs of an identical ligand–receptor pair cannot be attributed to the herb from which that ligand was catalogued.

Published experiments on β-sitosterol provide context but do not validate this study's proposed targets. Cao and colleagues investigated β-sitosterol in non-small-cell lung cancer models, not thyroid cancer.^17^ Chen and colleagues studied Prunella vulgaris components in inflammatory thyrocyte models relevant to Hashimoto's disease; these are not thyroid cancer models, and β-sitosterol showed no detectable affinity for HSP90AA1 in their reported binding assay.^18^ That result does not test CYP3A4 or CACNA1H, but it illustrates why network predictions require direct confirmation. Wnt and Keap1/Nrf2 signalling provide broader thyroid-cancer context,^19,20^ without establishing that these pathways were modulated in the present study.

Limitations include unresolved ingredient identity, incomplete compound–gene provenance, an unvalidated confidence-score threshold, and a broad disease set without subtype-specific or expression-based validation. Enrichment used default annotation backgrounds and did not establish target engagement or causal mechanisms. Database coverage, identifier mapping, shared ingredients, and term redundancy can influence the apparent profiles. Docking lacked validated receptor preparation, a fully documented search protocol, and suitable controls. No chemical standardisation, pharmacokinetic exposure, cellular response, or clinical outcome was measured. The immediate priorities are prescription reconciliation, human target-identity checks, reproducible enrichment with justified backgrounds, and controlled docking followed by target-specific experiments. The present evidence cannot support a recommendation for clinical use.

# Conclusion

The archived ingredient datasets contained genes annotated to thyroid-cancer-associated processes, but their overlap and enrichment profiles did not demonstrate therapeutic efficacy or superiority of XYT. Prescription inconsistency, a misidentified docking receptor, and unvalidated receptor preparation restrict the mechanistic conclusions. The results should be treated as provisional hypotheses pending data reconciliation and experimental verification.

# Conflict of Interest

None to declare, as stated in the supplied manuscript; all authors should reconfirm before submission.

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

The local review package contains recalculated set counts, the herb inventory, existing enrichment tables, figure source data, a docking audit, and a reproducible checking script with input-file hashes. [AUTHOR CONFIRM: public repository/accession or a justified data-access statement for submission, including the original compound–target source table and complete docking configurations.]

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
