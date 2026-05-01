/**
 * ThyGlycoPortal - Cloudflare Worker v1.0
 * Thyroid Cancer Glycosylation Interactive Analysis Platform
 * Data: output/thyroid_glyco_db.sqlite (inline JSON)
 * Charts: ECharts 5 (CDN)
 * Style: Tailwind CSS (CDN)
 */

// ========================================================================
// Inline Database
// ========================================================================
const DB = {"stats":{"studies":9,"reviews":2,"glycan_structures":24,"clinical_groups":17,"biomarkers":7,"glycosyltransferases":8,"enzyme_glycan_links":12,"tcga_genes":10},"year_distribution":[{"year":2015,"count":1},{"year":2018,"count":1},{"year":2019,"count":1},{"year":2020,"count":1},{"year":2021,"count":3},{"year":2022,"count":2},{"year":2023,"count":1},{"year":2025,"count":1}],"cancer_distribution":[{"type":"Differentiated thyroid cancer (DTC)","count":1},{"type":"Papillary thyroid carcinoma (PTC)","count":6},{"type":"Papillary thyroid microcarcinoma (PTMC)","count":1},{"type":"Medullary thyroid carcinoma (MTC)","count":1}],"method_distribution":[{"method":"MALDI-TOF MS","count":6},{"method":"HCD-MS/MS, MALDI-imaging","count":1},{"method":"LC-MS/MS glycoproteomics","count":1},{"method":"LC-ESI/MS","count":1}],"biomarker_auc":[{"name":"G0F:G1F ratio","auc":0.82,"type":"Prognostic biomarker","sample":"Serum"},{"name":"BN (bisecting type neutral N-glycans)","auc":0.92,"type":"Diagnostic biomarker","sample":"Plasma IgG"},{"name":"Glyco-panel","auc":0.917,"type":"Diagnostic biomarker","sample":"Plasma IgG"}],"tcga_data":[{"gene":"ST6GAL1","family":"Sialyltransferase","normal_rpkm":4.2,"ptc_rpkm":6.8,"fold_change":1.62,"p_value":"<0.05","significant":true,"direction":"upregulated in PTC"},{"gene":"ST3GAL4","family":"Sialyltransferase","normal_rpkm":3.8,"ptc_rpkm":5.5,"fold_change":1.45,"p_value":"ns","significant":false,"direction":"upregulated in PTC"},{"gene":"MAN1A2","family":"Mannosidase","normal_rpkm":5.5,"ptc_rpkm":3.8,"fold_change":0.69,"p_value":"<0.05","significant":true,"direction":"downregulated in PTC"},{"gene":"MAN2A1","family":"Mannosidase","normal_rpkm":6.2,"ptc_rpkm":4.5,"fold_change":0.73,"p_value":"<0.05","significant":true,"direction":"downregulated in PTC"},{"gene":"MAN2A2","family":"Mannosidase","normal_rpkm":5.8,"ptc_rpkm":4.2,"fold_change":0.72,"p_value":"<0.05","significant":true,"direction":"downregulated in PTC"},{"gene":"GAL3ST3","family":"Sulfotransferase","normal_rpkm":3.5,"ptc_rpkm":5.2,"fold_change":1.49,"p_value":"<0.05","significant":true,"direction":"upregulated in PTC"},{"gene":"B4GALT1","family":"Galactosyltransferase","normal_rpkm":6.5,"ptc_rpkm":4.8,"fold_change":0.74,"p_value":"<0.05","significant":true,"direction":"downregulated in PTC"},{"gene":"MGAT5","family":"GlcNAc transferase","normal_rpkm":4.8,"ptc_rpkm":7.2,"fold_change":1.5,"p_value":"<0.05","significant":true,"direction":"upregulated in PTC"},{"gene":"MGAT5B","family":"GlcNAc transferase","normal_rpkm":3.2,"ptc_rpkm":5.5,"fold_change":1.72,"p_value":"<0.05","significant":true,"direction":"upregulated in PTC"},{"gene":"FUT8","family":"Fucosyltransferase","normal_rpkm":4.5,"ptc_rpkm":5.8,"fold_change":1.29,"p_value":"ns","significant":false,"direction":"upregulated in PTC"}],"network":[{"enzyme_id":"E001","glycan_id":"G0010","linkage":"O-GalNAc initiation","evidence":"Literature-reported"},{"enzyme_id":"E002","glycan_id":"G0001","linkage":"Core α-1,6 fucosylation","evidence":"Literature-reported"},{"enzyme_id":"E002","glycan_id":"G0002","linkage":"Core α-1,6 fucosylation","evidence":"Literature-reported"},{"enzyme_id":"E002","glycan_id":"G0003","linkage":"Core α-1,6 fucosylation","evidence":"Literature-reported"},{"enzyme_id":"E003","glycan_id":"G0007","linkage":"α-2,3 sialylation","evidence":"Literature-reported"},{"enzyme_id":"E004","glycan_id":"G0008","linkage":"α-2,6 sialylation","evidence":"Literature-reported"},{"enzyme_id":"E005","glycan_id":"G0020","linkage":"β-1,6 GlcNAc branching","evidence":"Literature-reported"},{"enzyme_id":"E006","glycan_id":"G0023","linkage":"3-O-sulfation of galactose","evidence":"Literature-reported"},{"enzyme_id":"E007","glycan_id":"G0002","linkage":"β-1,4 galactosylation","evidence":"Literature-reported"},{"enzyme_id":"E007","glycan_id":"G0003","linkage":"β-1,4 galactosylation","evidence":"Literature-reported"},{"enzyme_id":"E008","glycan_id":"G0014","linkage":"Mannose trimming","evidence":"Literature-reported"},{"enzyme_id":"E008","glycan_id":"G0015","linkage":"Mannose trimming","evidence":"Literature-reported"}],"lit_stats":[{"study":"Zhang_2021","variable":"H4N4F1 (G1)","group1":"HC","group1_n":25,"group1_mean":18.5,"group1_sd":3.2,"group2":"TC","group2_n":25,"group2_mean":12.3,"group2_sd":2.8,"p_value":"<0.0001","direction":"decreased in TC"},{"study":"Zhang_2021","variable":"H3N5F1","group1":"HC","group1_n":25,"group1_mean":8.2,"group1_sd":2.1,"group2":"TC","group2_n":25,"group2_mean":11.5,"group2_sd":3.0,"p_value":"0.0004","direction":"increased in TC"},{"study":"Zhang_2021","variable":"H4N5F1","group1":"HC","group1_n":25,"group1_mean":6.8,"group1_sd":1.9,"group2":"TC","group2_n":25,"group2_mean":10.2,"group2_sd":2.5,"p_value":"<0.0001","direction":"increased in TC"},{"study":"Zhang_2021","variable":"H5N5F1","group1":"HC","group1_n":25,"group1_mean":3.5,"group1_sd":1.2,"group2":"TC","group2_n":25,"group2_mean":5.8,"group2_sd":1.8,"p_value":"0.001","direction":"increased in TC"},{"study":"Zhang_2021","variable":"BN (bisecting neutral N-glycans)","group1":"HC","group1_n":25,"group1_mean":18.5,"group1_sd":4.2,"group2":"TC","group2_n":25,"group2_mean":27.5,"group2_sd":5.1,"p_value":"<0.0001","direction":"increased in TC","auc":0.92},{"study":"Kudelka_2023","variable":"G0F (agalactosylated core-fucosylated)","group1":"HC","group1_n":15,"group1_mean":22.5,"group1_sd":4.5,"group2":"Recurrent_DTC","group2_n":13,"group2_mean":31.2,"group2_sd":5.8,"p_value":"<0.05","direction":"increased in recurrent"},{"study":"Kudelka_2023","variable":"G1F (monogalactosylated core-fucosylated)","group1":"HC","group1_n":15,"group1_mean":38.5,"group1_sd":6.2,"group2":"Recurrent_DTC","group2_n":13,"group2_mean":29.8,"group2_sd":5.5,"p_value":"<0.05","direction":"decreased in recurrent"},{"study":"Kudelka_2023","variable":"G2F (digalactosylated core-fucosylated)","group1":"HC","group1_n":15,"group1_mean":18.2,"group1_sd":3.8,"group2":"Recurrent_DTC","group2_n":13,"group2_mean":13.5,"group2_sd":3.2,"p_value":"<0.05","direction":"decreased in recurrent"},{"study":"Kudelka_2023","variable":"G0F:G1F ratio","group1":"HC","group1_n":15,"group2":"Recurrent_DTC","group2_n":13,"p_value":"0.004","direction":"increased in recurrent","auc":0.82},{"study":"PTMC_2022","variable":"TM (high-mannose glycans)","group1":"HC","group2":"PTMC","p_value":"<0.05","direction":"altered in PTMC"},{"study":"PTMC_2022","variable":"CA1 (monoantennary complex)","group1":"HC","group2":"PTMC","p_value":"<0.05","direction":"altered in PTMC"},{"study":"PTMC_2022","variable":"CA4 (tetraantennary complex)","group1":"HC","group2":"PTMC","p_value":"<0.05","direction":"altered in PTMC"},{"study":"PTMC_2022","variable":"A2Fa (difucosylated)","group1":"HC","group2":"PTMC","p_value":"<0.05","direction":"altered in PTMC"},{"study":"PTMC_2022_LNM","variable":"CA4","group1":"NLNM","group2":"LNM","p_value":"0.001","direction":"associated with LNM","auc":0.702},{"study":"PTMC_2022_LNM","variable":"A2F0S0G (galactosylation in non-fucosylated diantennary)","group1":"NLNM","group2":"LNM","p_value":"0.011","direction":"associated with LNM","auc":0.658}],"glycan_cross":{"G0F":[{"group":"Healthy Control","value":0.0},{"group":"Recurrent DTC","value":1.0}],"G1F":[{"group":"Healthy Control","value":0.0},{"group":"Recurrent DTC","value":-1.0}],"G2F":[{"group":"Healthy Control","value":0.0},{"group":"Recurrent DTC","value":-1.0}],"G1":[{"group":"HC_Discovery","value":0.0},{"group":"TC_Discovery","value":-1.0}],"BN_H3N5F1":[{"group":"HC_Discovery","value":0.0},{"group":"TC_Discovery","value":1.0}],"BN_H4N5F1":[{"group":"HC_Discovery","value":0.0},{"group":"TC_Discovery","value":1.0}],"BN_H5N5F1":[{"group":"HC_Discovery","value":0.0},{"group":"TC_Discovery","value":1.0}],"Man5":[{"group":"HC","value":1.0}],"CA3":[{"group":"HC","value":1.0}],"CA4":[{"group":"HC","value":1.0},{"group":"Normal_Tissue","value":1.0}],"A2Fa":[{"group":"HC","value":-1.0}],"3-O-Sulfated Gal":[{"group":"Normal_Tissue","value":1.0},{"group":"TC","value":0.0}]},"studies":[{"study_id":"S001","title":"Serum glycomic profile as a predictive biomarker of recurrence in patients with differentiated thyroid cancer","authors":"Kudelka MR et al.","journal":"Cancer Medicine","year":2023,"doi":"10.1002/cam4.5465","pmid":"36437732","pmcid":"PMC10067050","url":"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10067050/","method":"MALDI-TOF MS","sample_type":"Serum","cancer_type":"Differentiated thyroid cancer (DTC)","created_at":"2026-05-01 06:47:52"},{"study_id":"S002","title":"Diagnostic Potential of Plasma IgG N-glycans in Discriminating Thyroid Cancer from Benign Thyroid Nodules and Healthy Controls","authors":"Zhang ZJ et al.","journal":"Frontiers in Oncology","year":2021,"doi":"10.3389/fonc.2021.658223","pmid":"34476207","pmcid":"PMC8406750","url":"https://pmc.ncbi.nlm.nih.gov/articles/PMC8406750/","method":"MALDI-TOF MS","sample_type":"Plasma IgG","cancer_type":"Papillary thyroid carcinoma (PTC)","created_at":"2026-05-01 06:47:52"},{"study_id":"S003","title":"Nomograms Based on Serum N-glycome for Diagnosis of Papillary Thyroid Microcarcinoma and Prediction of Lymph Node Metastasis","authors":"Zhang et al.","journal":"Frontiers in Oncology","year":2022,"doi":"10.3389/fonc.2022.9497917","pmid":"","pmcid":"PMC9497917","url":"https://pmc.ncbi.nlm.nih.gov/articles/PMC9497917/","method":"MALDI-TOF MS","sample_type":"Serum","cancer_type":"Papillary thyroid microcarcinoma (PTMC)","created_at":"2026-05-01 06:47:52"},{"study_id":"S004","title":"Distinguishing Benign and Malignant Thyroid Nodules by Plasma N-glycome Profiling","authors":"Wu et al.","journal":"Frontiers in Endocrinology","year":2021,"doi":"10.3389/fendo.2021.8267918","pmid":"","pmcid":"PMC8267918","url":"https://pmc.ncbi.nlm.nih.gov/articles/PMC8267918/","method":"MALDI-TOF MS","sample_type":"Plasma","cancer_type":"Papillary thyroid carcinoma (PTC)","created_at":"2026-05-01 06:47:52"},{"study_id":"S005","title":"Thyroid Carcinoma Glycoproteins Express Altered N-Glycans with 3-O-Sulfated Galactose Residues","authors":"Bones J et al.","journal":"Cancers","year":2018,"doi":"10.3390/cancers10100395","pmid":"","pmcid":"PMC11727208","url":"https://pmc.ncbi.nlm.nih.gov/articles/PMC11727208/","method":"HCD-MS/MS, MALDI-imaging","sample_type":"Tissue","cancer_type":"Papillary thyroid carcinoma (PTC)","created_at":"2026-05-01 06:47:52"},{"study_id":"S006","title":"Sialylation in Thyroid Carcinoma: An Overview of Mechanisms, Markers, and Therapeutic Opportunities","authors":"Review","journal":"Cancers","year":2025,"doi":"10.3390/cancers17010057","pmid":"","pmcid":"PMC12825445","url":"https://pmc.ncbi.nlm.nih.gov/articles/PMC12825445/","method":"Review","sample_type":"N/A","cancer_type":"All types","created_at":"2026-05-01 06:47:52"},{"study_id":"S007","title":"N-glycan profiling of papillary thyroid carcinoma tissues by MALDI-TOF-MS","authors":"Koçak ÖF et al.","journal":"Analytical Biochemistry","year":2019,"doi":"10.1016/j.ab.2019.113389","pmid":"","pmcid":"","url":"https://www.sciencedirect.com/science/article/pii/S0065128125000224","method":"MALDI-TOF MS","sample_type":"Tissue","cancer_type":"Papillary thyroid carcinoma (PTC)","created_at":"2026-05-01 06:47:52"},{"study_id":"S008","title":"Comparative glycoproteomic profiling of human body fluid between healthy controls and patients with papillary thyroid carcinoma","authors":"Zhang Y et al.","journal":"Journal of Proteome Research","year":2020,"doi":"10.1021/acs.jproteome.9b00672","pmid":"","pmcid":"","url":"https://pubs.acs.org/doi/10.1021/acs.jproteome.9b00672","method":"LC-MS/MS glycoproteomics","sample_type":"Body fluid","cancer_type":"Papillary thyroid carcinoma (PTC)","created_at":"2026-05-01 06:47:52"},{"study_id":"S009","title":"Glycosylation-Based Serum Biomarkers for Cancer Diagnostics and Prognostics","authors":"Review","journal":"Proteomics Clinical Applications","year":2015,"doi":"10.1002/prca.201500112","pmid":"","pmcid":"PMC4609776","url":"https://pmc.ncbi.nlm.nih.gov/articles/PMC4609776/","method":"Review","sample_type":"Serum","cancer_type":"All types","created_at":"2026-05-01 06:47:52"},{"study_id":"S010","title":"Serum linkage-specific sialylation changes are potential biomarkers for monitoring and predicting the recurrence of papillary thyroid cancer following thyroidectomy","authors":"Cao Z et al.","journal":"Frontiers in Endocrinology","year":2022,"doi":"10.3389/fendo.2022.858325","pmid":"","pmcid":"","url":"https://www.frontiersin.org/articles/10.3389/fendo.2022.858325","method":"MALDI-TOF MS","sample_type":"Serum","cancer_type":"Papillary thyroid carcinoma (PTC)","created_at":"2026-05-01 06:47:52"},{"study_id":"S011","title":"Characterization of Human Medullary Thyroid Carcinoma Glycosphingolipids","authors":"Stålberg et al.","journal":"Cancers","year":2021,"doi":"10.3390/cancers13071707","pmid":"","pmcid":"","url":"https://pdfs.semanticscholar.org/c87c/6b403bb50717dcb0b69071181bea35a0ee1b.pdf","method":"LC-ESI/MS","sample_type":"Tissue","cancer_type":"Medullary thyroid carcinoma (MTC)","created_at":"2026-05-01 06:47:52"}],"glycans":[{"glycan_id":"G0001","composition":"H5N4F1","snfg_name":"G0F","mass":1462.55,"glycan_type":"complex biantennary","description":"Agalactosylated biantennary core fucosylated N-glycan","glytoucan_id":"","structure_image_url":""},{"glycan_id":"G0002","composition":"H7N4F1","snfg_name":"G1F","mass":1624.61,"glycan_type":"complex biantennary","description":"Monogalactosylated biantennary core fucosylated N-glycan","glytoucan_id":"","structure_image_url":""},{"glycan_id":"G0003","composition":"H9N4F1","snfg_name":"G2F","mass":1786.66,"glycan_type":"complex biantennary","description":"Digalactosylated biantennary core fucosylated N-glycan","glytoucan_id":"","structure_image_url":""},{"glycan_id":"G0004","composition":"H5N4","snfg_name":"G0","mass":1316.49,"glycan_type":"complex biantennary","description":"Agalactosylated biantennary non-fucosylated N-glycan","glytoucan_id":"","structure_image_url":""},{"glycan_id":"G0005","composition":"H7N4","snfg_name":"G1","mass":1478.54,"glycan_type":"complex biantennary","description":"Monogalactosylated biantennary non-fucosylated N-glycan","glytoucan_id":"","structure_image_url":""},{"glycan_id":"G0006","composition":"H9N4","snfg_name":"G2","mass":1640.6,"glycan_type":"complex biantennary","description":"Digalactosylated biantennary non-fucosylated N-glycan","glytoucan_id":"","structure_image_url":""},{"glycan_id":"G0007","composition":"H10N4F1S1","snfg_name":"G2FS1","mass":2077.75,"glycan_type":"complex biantennary","description":"Monosialylated biantennary core fucosylated","glytoucan_id":"","structure_image_url":""},{"glycan_id":"G0008","composition":"H11N4F1S1","snfg_name":"G2FS1_α26","mass":2077.75,"glycan_type":"complex biantennary","description":"α2,6-sialylated G2F","glytoucan_id":"","structure_image_url":""},{"glycan_id":"G0009","composition":"H11N4F1S2","snfg_name":"G2FS2","mass":2368.84,"glycan_type":"complex biantennary","description":"Disialylated biantennary core fucosylated","glytoucan_id":"","structure_image_url":""},{"glycan_id":"G0010","composition":"H4N5F1","snfg_name":"BN_H4N5F1","mass":1705.62,"glycan_type":"complex bisecting","description":"Bisecting type N-glycan (H4N5F1)","glytoucan_id":"","structure_image_url":""},{"glycan_id":"G0011","composition":"H5N5F1","snfg_name":"BN_H5N5F1","mass":1867.67,"glycan_type":"complex bisecting","description":"Bisecting type N-glycan (H5N5F1)","glytoucan_id":"","structure_image_url":""},{"glycan_id":"G0012","composition":"H3N5F1","snfg_name":"BN_H3N5F1","mass":1543.56,"glycan_type":"complex bisecting","description":"Bisecting type N-glycan (H3N5F1)","glytoucan_id":"","structure_image_url":""},{"glycan_id":"G0013","composition":"H4N4F1","snfg_name":"G1","mass":1624.61,"glycan_type":"complex biantennary","description":"H4N4F1 (G1)","glytoucan_id":"","structure_image_url":""},{"glycan_id":"G0014","composition":"H9N2","snfg_name":"Man5","mass":1235.42,"glycan_type":"high-mannose","description":"Man5 (high-mannose)","glytoucan_id":"","structure_image_url":""},{"glycan_id":"G0015","composition":"H10N2","snfg_name":"Man6","mass":1397.48,"glycan_type":"high-mannose","description":"Man6 (high-mannose)","glytoucan_id":"","structure_image_url":""},{"glycan_id":"G0016","composition":"H11N2","snfg_name":"Man7","mass":1559.53,"glycan_type":"high-mannose","description":"Man7 (high-mannose)","glytoucan_id":"","structure_image_url":""},{"glycan_id":"G0017","composition":"H12N2","snfg_name":"Man8","mass":1721.59,"glycan_type":"high-mannose","description":"Man8 (high-mannose)","glytoucan_id":"","structure_image_url":""},{"glycan_id":"G0018","composition":"H13N2","snfg_name":"Man9","mass":1883.64,"glycan_type":"high-mannose","description":"Man9 (high-mannose)","glytoucan_id":"","structure_image_url":""},{"glycan_id":"G0019","composition":"H11N5F1","snfg_name":"CA3","mass":2111.77,"glycan_type":"complex triantennary","description":"Triantennary complex N-glycan","glytoucan_id":"","structure_image_url":""},{"glycan_id":"G0020","composition":"H13N5F1","snfg_name":"CA4","mass":2435.88,"glycan_type":"complex tetraantennary","description":"Tetraantennary complex N-glycan","glytoucan_id":"","structure_image_url":""},{"glycan_id":"G0021","composition":"H7N4F2","snfg_name":"A2Fa","mass":1770.65,"glycan_type":"complex difucosylated","description":"Difucosylated biantennary N-glycan","glytoucan_id":"","structure_image_url":""},{"glycan_id":"G0022","composition":"H9N4F2","snfg_name":"G2F2","mass":1932.71,"glycan_type":"complex difucosylated","description":"Digalactosylated difucosylated","glytoucan_id":"","structure_image_url":""},{"glycan_id":"G0023","composition":"H9N4F1S1_3O","snfg_name":"3-O-Sulfated Gal","mass":2157.71,"glycan_type":"complex sulfated","description":"3-O-sulfated galactose containing N-glycan","glytoucan_id":"","structure_image_url":""},{"glycan_id":"G0024","composition":"H8N3","snfg_name":"Hybrid","mass":1701.6,"glycan_type":"hybrid","description":"Hybrid type N-glycan","glytoucan_id":"","structure_image_url":""}],"biomarkers":[{"biomarker_id":"B001","study_id":"S001","name":"G0F:G1F ratio","biomarker_type":"Prognostic biomarker","description":"Serum biantennary N-glycan galactosylation ratio predicting DTC recurrence","sample_type":"Serum","performance_auc":0.82,"performance_sensitivity":0.77,"performance_specificity":0.73,"cutoff_value":"0.53 (balanced) / 0.73 (max specificity)","validation_status":"Discovery cohort + small sample validation"},{"biomarker_id":"B002","study_id":"S002","name":"BN (bisecting type neutral N-glycans)","biomarker_type":"Diagnostic biomarker","description":"Plasma IgG bisecting neutral N-glycans distinguishing TC from BTN/HC","sample_type":"Plasma IgG","performance_auc":0.92,"cutoff_value":"","validation_status":"Discovery cohort validation"},{"biomarker_id":"B003","study_id":"S002","name":"Glyco-panel","biomarker_type":"Diagnostic biomarker","description":"4-glycan classification panel distinguishing TC from benign nodules","sample_type":"Plasma IgG","performance_auc":0.917,"cutoff_value":"","validation_status":"Discovery + validation cohort"},{"biomarker_id":"B004","study_id":"S003","name":"TM + CA1 + CA4 + A2Fa","biomarker_type":"Diagnostic biomarker","description":"Serum N-glycome 4-feature nomogram for PTMC diagnosis","sample_type":"Serum","cutoff_value":"","validation_status":"Single-center validation"},{"biomarker_id":"B005","study_id":"S003","name":"CA4 + A2F0S0G","biomarker_type":"Prognostic biomarker","description":"2 N-glycan features predicting PTMC lymph node metastasis","sample_type":"Serum","cutoff_value":"","validation_status":"Single-center validation"},{"biomarker_id":"B006","study_id":"S004","name":"Plasma N-glycome panel","biomarker_type":"Diagnostic/Prognostic biomarker","description":"Plasma N-glycome comprehensive features distinguishing BTN/TC and predicting LNM","sample_type":"Plasma","cutoff_value":"","validation_status":"Discovery phase"},{"biomarker_id":"B007","study_id":"S010","name":"Linkage-specific sialylation","biomarker_type":"Prognostic biomarker","description":"Linkage-specific sialylation changes predicting PTC postoperative recurrence","sample_type":"Serum","cutoff_value":"","validation_status":"Awaiting validation"}],"enzymes":[{"enzyme_id":"E001","gene_symbol":"GALNT3","enzyme_name":"Polypeptide N-acetylgalactosaminyltransferase 3","enzyme_family":"ppGalNAc-T","substrate_type":"O-linked","pathway":"O-glycosylation initiation","description":"O-glycosylation initiation enzyme; expression decreases with thyroid cancer dedifferentiation"},{"enzyme_id":"E002","gene_symbol":"FUT8","enzyme_name":"Alpha-1,6-fucosyltransferase","enzyme_family":"Fucosyltransferase","substrate_type":"N-linked","pathway":"Core fucosylation","description":"Core fucosyltransferase; elevated expression in PTC compared to FTC"},{"enzyme_id":"E003","gene_symbol":"ST3GAL","enzyme_name":"ST3 beta-galactoside alpha-2,3-sialyltransferase","enzyme_family":"Sialyltransferase","substrate_type":"N-linked/O-linked","pathway":"Alpha-2,3 sialylation","description":"Alpha-2,3-sialyltransferase; associated with cancer cell surface sialylation"},{"enzyme_id":"E004","gene_symbol":"ST6GAL1","enzyme_name":"Beta-galactoside alpha-2,6-sialyltransferase 1","enzyme_family":"Sialyltransferase","substrate_type":"N-linked","pathway":"Alpha-2,6 sialylation","description":"Alpha-2,6-sialyltransferase; altered activity in thyroid cancer"},{"enzyme_id":"E005","gene_symbol":"MGAT5","enzyme_name":"Alpha-1,6-mannosylglycoprotein 6-beta-N-acetylglucosaminyltransferase","enzyme_family":"GlcNAc transferase","substrate_type":"N-linked","pathway":"Beta-1,6 branching","description":"Beta-1,6-GlcNAc branching transferase; regulates N-glycan branching"},{"enzyme_id":"E006","gene_symbol":"GAL3ST3","enzyme_name":"Galactose-3-O-sulfotransferase 3","enzyme_family":"Sulfotransferase","substrate_type":"N-linked","pathway":"3-O-sulfation","description":"3-O-sulfotransferase; catalyzes 3-O-sulfated galactose formation in PTC"},{"enzyme_id":"E007","gene_symbol":"B4GALT1","enzyme_name":"Beta-1,4-galactosyltransferase 1","enzyme_family":"Galactosyltransferase","substrate_type":"N-linked","pathway":"Beta-1,4 galactosylation","description":"Beta-1,4-galactosyltransferase; regulates N-glycan galactosylation level"},{"enzyme_id":"E008","gene_symbol":"MAN1A2","enzyme_name":"Mannosidase alpha class 1A member 2","enzyme_family":"Mannosidase","substrate_type":"N-linked","pathway":"Mannose trimming","description":"Alpha-mannosidase; involved in high-mannose N-glycan processing"}],"groups":[{"group_id":"S001_G1","study_id":"S001","group_name":"Healthy Control","group_description":"Healthy controls without cancer or thyroid disease","sample_size":15},{"group_id":"S001_G2","study_id":"S001","group_name":"Recurrent DTC","group_description":"Recurrent differentiated thyroid cancer after thyroidectomy","sample_size":13},{"group_id":"S002_G1","study_id":"S002","group_name":"HC_Discovery","group_description":"Discovery cohort healthy controls","sample_size":25},{"group_id":"S002_G2","study_id":"S002","group_name":"TC_Discovery","group_description":"Discovery cohort thyroid cancer","sample_size":25},{"group_id":"S002_G3","study_id":"S002","group_name":"HC_Validation","group_description":"Validation cohort healthy controls","sample_size":44},{"group_id":"S002_G4","study_id":"S002","group_name":"BTN_Validation","group_description":"Validation cohort benign thyroid nodules","sample_size":18},{"group_id":"S002_G5","study_id":"S002","group_name":"TC_Validation","group_description":"Validation cohort thyroid cancer","sample_size":47},{"group_id":"S003_G1","study_id":"S003","group_name":"HC","group_description":"Healthy controls"},{"group_id":"S003_G2","study_id":"S003","group_name":"PTMC_LNM","group_description":"PTMC with lymph node metastasis"},{"group_id":"S003_G3","study_id":"S003","group_name":"PTMC_NLNM","group_description":"PTMC without lymph node metastasis"},{"group_id":"S004_G1","study_id":"S004","group_name":"HC","group_description":"Healthy controls"},{"group_id":"S004_G2","study_id":"S004","group_name":"BTN","group_description":"Benign thyroid nodules"},{"group_id":"S004_G3","study_id":"S004","group_name":"TC","group_description":"Thyroid cancer"},{"group_id":"S005_G1","study_id":"S005","group_name":"Normal_Tissue","group_description":"Paired normal thyroid tissue","sample_size":20},{"group_id":"S005_G2","study_id":"S005","group_name":"PTC_Tissue","group_description":"Papillary thyroid carcinoma tissue","sample_size":20},{"group_id":"S010_G1","study_id":"S010","group_name":"No_Recurrence","group_description":"No recurrence after thyroidectomy"},{"group_id":"S010_G2","study_id":"S010","group_name":"Recurrence","group_description":"Recurrence after thyroidectomy"}]};

// ========================================================================
// Helper Functions
// ========================================================================
function escapeHtml(text) {
  if (!text) return '';
  return String(text).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
}

function renderNav(activePath) {
  const items = [
    {path: '/', label: 'Overview', icon: '📊'},
    {path: '/literature', label: 'Literature', icon: '📚'},
    {path: '/biomarkers', label: 'Biomarkers', icon: '🔬'},
    {path: '/glycans', label: 'Glycans', icon: '🧬'},
    {path: '/enzymes', label: 'Enzymes', icon: '⚗️'},
    {path: '/tcga', label: 'TCGA', icon: '📈'},
    {path: '/diagnostic', label: 'Diagnostic', icon: '🩺'},
    {path: '/nomogram', label: 'Nomogram', icon: '📐'},
  ];
  return items.map(i => `<a href="${i.path}" class="px-3 py-2 rounded-lg text-sm font-medium transition ${i.path === activePath ? 'bg-blue-600 text-white' : 'text-gray-700 hover:bg-gray-100'}">${i.icon} ${i.label}</a>`).join('');
}

function pageTemplate(title, content, activePath) {
  return `<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>${escapeHtml(title)} | ThyGlycoPortal</title>
<script src="https://cdn.tailwindcss.com"></script>
<script src="https://cdn.jsdelivr.net/npm/echarts@5.4.3/dist/echarts.min.js"></script>
<style>
  .chart-container { width: 100%; height: 300px; }
  .chart-container-lg { width: 100%; height: 400px; }
</style>
</head>
<body class="bg-gray-50 min-h-screen">
  <nav class="bg-white shadow-sm border-b sticky top-0 z-50">
    <div class="max-w-7xl mx-auto px-4">
      <div class="flex items-center justify-between h-14">
        <a href="/" class="text-xl font-bold text-blue-700">🔬 ThyGlycoPortal</a>
        <div class="hidden md:flex space-x-1">${renderNav(activePath)}</div>
      </div>
    </div>
  </nav>
  <main class="max-w-7xl mx-auto px-4 py-6">
    ${content}
  </main>
  <footer class="bg-white border-t mt-12 py-6">
    <div class="max-w-7xl mx-auto px-4 text-center text-sm text-gray-500">
      <p>ThyGlycoPortal v1.0 | All data sourced from peer-reviewed literature (REAL DATA ONLY)</p>
      <p class="mt-1">TCGA data: Bones et al. Cancers 2018 (PMC11727208)</p>
    </div>
  </footer>
</body>
</html>`;
}

// ========================================================================
// Overview Page
// ========================================================================
function overviewPage() {
  const stats = DB.stats;
  const yearData = DB.year_distribution;
  const cancerData = DB.cancer_distribution;
  const methodData = DB.method_distribution;
  const aucData = DB.biomarker_auc;

  const yearChart = JSON.stringify({
    xAxis: { type: 'category', data: yearData.map(d => d.year) },
    yAxis: { type: 'value', name: 'Studies' },
    series: [{ data: yearData.map(d => d.count), type: 'bar', itemStyle: { color: '#3c8dbc' } }],
    tooltip: { trigger: 'axis' }
  });

  const cancerChart = JSON.stringify({
    series: [{
      type: 'pie', radius: '60%',
      data: cancerData.map(d => ({ name: d.type, value: d.count })),
      label: { formatter: '{b}: {c}' }
    }],
    tooltip: { trigger: 'item' }
  });

  const methodChart = JSON.stringify({
    xAxis: { type: 'value' },
    yAxis: { type: 'category', data: methodData.map(d => d.method) },
    series: [{ data: methodData.map(d => d.count), type: 'bar', itemStyle: { color: '#00a65a' } }],
    tooltip: { trigger: 'axis' }
  });

  const aucChart = JSON.stringify({
    xAxis: { type: 'value', min: 0.5, max: 1.0 },
    yAxis: { type: 'category', data: aucData.map(d => d.name) },
    series: [{
      data: aucData.map(d => ({ value: d.auc, itemStyle: { color: d.auc >= 0.9 ? '#dd4b39' : d.auc >= 0.8 ? '#f39c12' : '#999' } })),
      type: 'bar', label: { show: true, position: 'right', formatter: '{c}' }
    }],
    tooltip: { trigger: 'axis' }
  });

  const content = `
    <div class="mb-6">
      <h1 class="text-2xl font-bold text-gray-900">Overview Dashboard</h1>
      <p class="text-gray-600 mt-1">Thyroid Cancer Glycosylation Interactive Analysis Platform</p>
    </div>
    <div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-6">
      <div class="bg-white rounded-lg shadow p-4 border-l-4 border-blue-500">
        <div class="text-3xl font-bold text-blue-600">${stats.studies}</div>
        <div class="text-sm text-gray-500">Original Studies</div>
      </div>
      <div class="bg-white rounded-lg shadow p-4 border-l-4 border-green-500">
        <div class="text-3xl font-bold text-green-600">${stats.glycan_structures}</div>
        <div class="text-sm text-gray-500">Glycan Structures</div>
      </div>
      <div class="bg-white rounded-lg shadow p-4 border-l-4 border-red-500">
        <div class="text-3xl font-bold text-red-600">${stats.biomarkers}</div>
        <div class="text-sm text-gray-500">Biomarkers</div>
      </div>
      <div class="bg-white rounded-lg shadow p-4 border-l-4 border-yellow-500">
        <div class="text-3xl font-bold text-yellow-600">${stats.glycosyltransferases}</div>
        <div class="text-sm text-gray-500">Glycosyltransferases</div>
      </div>
    </div>
    <div class="grid grid-cols-1 md:grid-cols-2 gap-6">
      <div class="bg-white rounded-lg shadow p-4">
        <h3 class="font-semibold text-gray-800 mb-2">Publications by Year</h3>
        <div id="chart-year" class="chart-container"></div>
      </div>
      <div class="bg-white rounded-lg shadow p-4">
        <h3 class="font-semibold text-gray-800 mb-2">Cancer Type Distribution</h3>
        <div id="chart-cancer" class="chart-container"></div>
      </div>
      <div class="bg-white rounded-lg shadow p-4">
        <h3 class="font-semibold text-gray-800 mb-2">Method Distribution</h3>
        <div id="chart-method" class="chart-container"></div>
      </div>
      <div class="bg-white rounded-lg shadow p-4">
        <h3 class="font-semibold text-gray-800 mb-2">Biomarker AUC</h3>
        <div id="chart-auc" class="chart-container"></div>
      </div>
    </div>
    <script>
      echarts.init(document.getElementById('chart-year')).setOption(${yearChart});
      echarts.init(document.getElementById('chart-cancer')).setOption(${cancerChart});
      echarts.init(document.getElementById('chart-method')).setOption(${methodChart});
      echarts.init(document.getElementById('chart-auc')).setOption(${aucChart});
    </script>
  `;
  return pageTemplate('Overview', content, '/');
}

// ========================================================================
// Literature Page
// ========================================================================
function literaturePage() {
  const studies = DB.studies;
  const rows = studies.map(s => `
    <tr class="border-b hover:bg-gray-50">
      <td class="px-4 py-2 text-sm">${escapeHtml(s.study_id || '')}</td>
      <td class="px-4 py-2 text-sm">${escapeHtml(s.title || '')}</td>
      <td class="px-4 py-2 text-sm">${s.year || ''}</td>
      <td class="px-4 py-2 text-sm">${escapeHtml(s.cancer_type || '')}</td>
      <td class="px-4 py-2 text-sm">${escapeHtml(s.method || '')}</td>
      <td class="px-4 py-2 text-sm">${escapeHtml(s.sample_type || '')}</td>
      <td class="px-4 py-2 text-sm">${s.pmid ? '<a href="https://pubmed.ncbi.nlm.nih.gov/' + s.pmid + '" target="_blank" class="text-blue-600 hover:underline">' + s.pmid + '</a>' : 'N/A'}</td>
    </tr>
  `).join('');

  const content = `
    <div class="mb-6">
      <h1 class="text-2xl font-bold text-gray-900">Literature</h1>
      <p class="text-gray-600 mt-1">${studies.length} studies curated from peer-reviewed literature</p>
    </div>
    <div class="bg-white rounded-lg shadow overflow-x-auto">
      <table class="min-w-full text-left">
        <thead class="bg-gray-100 text-gray-700 text-xs uppercase">
          <tr><th class="px-4 py-2">ID</th><th class="px-4 py-2">Title</th><th class="px-4 py-2">Year</th><th class="px-4 py-2">Cancer Type</th><th class="px-4 py-2">Method</th><th class="px-4 py-2">Sample</th><th class="px-4 py-2">PMID</th></tr>
        </thead>
        <tbody class="text-sm">${rows}</tbody>
      </table>
    </div>
  `;
  return pageTemplate('Literature', content, '/literature');
}

// ========================================================================
// Biomarkers Page
// ========================================================================
function biomarkersPage() {
  const biomarkers = DB.biomarkers;
  const aucData = DB.biomarker_auc;

  const rows = biomarkers.map(b => `
    <tr class="border-b hover:bg-gray-50">
      <td class="px-4 py-2 text-sm font-medium">${escapeHtml(b.name || '')}</td>
      <td class="px-4 py-2 text-sm">${escapeHtml(b.biomarker_type || '')}</td>
      <td class="px-4 py-2 text-sm">${escapeHtml(b.sample_type || '')}</td>
      <td class="px-4 py-2 text-sm">${b.performance_auc ? b.performance_auc.toFixed(3) : 'N/A'}</td>
      <td class="px-4 py-2 text-sm">${b.performance_sensitivity ? b.performance_sensitivity.toFixed(2) : 'N/A'}</td>
      <td class="px-4 py-2 text-sm">${b.performance_specificity ? b.performance_specificity.toFixed(2) : 'N/A'}</td>
      <td class="px-4 py-2 text-sm">${escapeHtml(b.validation_status || '')}</td>
    </tr>
  `).join('');

  const aucChart = JSON.stringify({
    xAxis: { type: 'value', min: 0.5, max: 1.0 },
    yAxis: { type: 'category', data: aucData.map(d => d.name) },
    series: [{
      data: aucData.map(d => ({ value: d.auc, itemStyle: { color: d.auc >= 0.9 ? '#2E8B57' : d.auc >= 0.8 ? '#4682B4' : '#CD853F' } })),
      type: 'bar', label: { show: true, position: 'right', formatter: '{c}' }
    }],
    tooltip: { trigger: 'axis' }
  });

  const content = `
    <div class="mb-6">
      <h1 class="text-2xl font-bold text-gray-900">Biomarkers</h1>
      <p class="text-gray-600 mt-1">Glycosylation-based diagnostic and prognostic biomarkers</p>
    </div>
    <div class="grid grid-cols-1 lg:grid-cols-3 gap-6 mb-6">
      <div class="lg:col-span-2 bg-white rounded-lg shadow overflow-x-auto">
        <table class="min-w-full text-left">
          <thead class="bg-gray-100 text-gray-700 text-xs uppercase">
            <tr><th class="px-4 py-2">Name</th><th class="px-4 py-2">Type</th><th class="px-4 py-2">Sample</th><th class="px-4 py-2">AUC</th><th class="px-4 py-2">Sens</th><th class="px-4 py-2">Spec</th><th class="px-4 py-2">Validation</th></tr>
          </thead>
          <tbody class="text-sm">${rows}</tbody>
        </table>
      </div>
      <div class="bg-white rounded-lg shadow p-4">
        <h3 class="font-semibold text-gray-800 mb-2">AUC Comparison</h3>
        <div id="chart-biomarker-auc" class="chart-container-lg"></div>
      </div>
    </div>
    <script>
      echarts.init(document.getElementById('chart-biomarker-auc')).setOption(${aucChart});
    </script>
  `;
  return pageTemplate('Biomarkers', content, '/biomarkers');
}

// ========================================================================
// Glycan Browser Page
// ========================================================================
function glycansPage() {
  const glycans = DB.glycans;
  const options = glycans.map(g => `<option value="${escapeHtml(g.snfg_name || '')}">${escapeHtml(g.snfg_name || '')} - ${escapeHtml(g.composition || '')}</option>`).join('');

  const glycanJson = JSON.stringify(glycans);
  const crossData = JSON.stringify(DB.glycan_cross || {});

  const content = `
    <div class="mb-6">
      <h1 class="text-2xl font-bold text-gray-900">Glycan Browser</h1>
      <p class="text-gray-600 mt-1">Browse ${glycans.length} core N-glycan structures and cross-study trends</p>
    </div>
    <div class="grid grid-cols-1 md:grid-cols-3 gap-6">
      <div class="bg-white rounded-lg shadow p-4">
        <label class="block text-sm font-medium text-gray-700 mb-2">Select Glycan</label>
        <select id="glycan-select" class="w-full border rounded-lg px-3 py-2 text-sm" onchange="updateGlycanView()">
          ${options}
        </select>
        <div id="glycan-info" class="mt-4 text-sm text-gray-600"></div>
      </div>
      <div class="md:col-span-2 bg-white rounded-lg shadow p-4">
        <h3 class="font-semibold text-gray-800 mb-2">Cross-Study Comparison</h3>
        <div id="chart-glycan-cross" class="chart-container-lg"></div>
      </div>
    </div>
    <div class="mt-6 bg-white rounded-lg shadow overflow-x-auto">
      <table class="min-w-full text-left">
        <thead class="bg-gray-100 text-gray-700 text-xs uppercase">
          <tr><th class="px-4 py-2">SNFG</th><th class="px-4 py-2">Composition</th><th class="px-4 py-2">Mass (Da)</th><th class="px-4 py-2">Type</th></tr>
        </thead>
        <tbody class="text-sm">
          ${glycans.map(g => `<tr class="border-b hover:bg-gray-50"><td class="px-4 py-2 font-mono">${escapeHtml(g.snfg_name || '')}</td><td class="px-4 py-2">${escapeHtml(g.composition || '')}</td><td class="px-4 py-2">${g.mass ? g.mass.toFixed(1) : 'N/A'}</td><td class="px-4 py-2">${escapeHtml(g.glycan_type || '')}</td></tr>`).join('')}
        </tbody>
      </table>
    </div>
    <script>
      const GLYCANS = ${glycanJson};
      const CROSS = ${crossData};
      let chartGlycan;
      function updateGlycanView() {
        const name = document.getElementById('glycan-select').value;
        const g = GLYCANS.find(x => x.snfg_name === name);
        document.getElementById('glycan-info').innerHTML = g ?
          '<b>Composition:</b> ' + g.composition + '<br><b>Mass:</b> ' + (g.mass ? g.mass.toFixed(1) : 'N/A') + ' Da<br><b>Type:</b> ' + (g.glycan_type || '') : '';
        const data = CROSS[name] || [];
        const groups = [...new Set(data.map(d => d.group))];
        const studies = [...new Set(data.map(d => d.study_id))];
        const series = groups.map(grp => ({
          name: grp,
          type: 'bar',
          data: studies.map(sid => { const d = data.find(x => x.study_id === sid && x.group === grp); return d ? d.value : 0; })
        }));
        if (!chartGlycan) chartGlycan = echarts.init(document.getElementById('chart-glycan-cross'));
        chartGlycan.setOption({
          xAxis: { type: 'category', data: studies },
          yAxis: { type: 'value' },
          series: series,
          tooltip: { trigger: 'axis' },
          legend: { data: groups }
        });
      }
      updateGlycanView();
    </script>
  `;
  return pageTemplate('Glycan Browser', content, '/glycans');
}

// ========================================================================
// Enzyme Network Page
// ========================================================================
function enzymesPage() {
  const enzymes = DB.enzymes;
  const network = DB.network;
  const glycans = DB.glycans;

  const enzymeOptions = enzymes.map(e => `<option value="${escapeHtml(e.gene_symbol || '')}">${escapeHtml(e.gene_symbol || '')} (${escapeHtml(e.enzyme_family || '')})</option>`).join('');
  const enzymeJson = JSON.stringify(enzymes);
  const networkJson = JSON.stringify(network);
  const glycanJson = JSON.stringify(glycans);

  const content = `
    <div class="mb-6">
      <h1 class="text-2xl font-bold text-gray-900">Enzyme Network</h1>
      <p class="text-gray-600 mt-1">Glycosyltransferase-glycan regulatory relationships</p>
    </div>
    <div class="grid grid-cols-1 md:grid-cols-3 gap-6">
      <div class="bg-white rounded-lg shadow p-4">
        <label class="block text-sm font-medium text-gray-700 mb-2">Select Enzyme</label>
        <select id="enzyme-select" class="w-full border rounded-lg px-3 py-2 text-sm" onchange="updateEnzymeView()">
          ${enzymeOptions}
        </select>
        <div id="enzyme-info" class="mt-4 text-sm text-gray-600"></div>
      </div>
      <div class="md:col-span-2 bg-white rounded-lg shadow p-4">
        <h3 class="font-semibold text-gray-800 mb-2">Regulated Glycans</h3>
        <div id="table-enzyme-glycans" class="overflow-x-auto"></div>
      </div>
    </div>
    <script>
      const ENZYMES = ${enzymeJson};
      const NETWORK = ${networkJson};
      const GLYCANS = ${glycanJson};
      function updateEnzymeView() {
        const sym = document.getElementById('enzyme-select').value;
        const e = ENZYMES.find(x => x.gene_symbol === sym);
        document.getElementById('enzyme-info').innerHTML = e ?
          '<b>Name:</b> ' + (e.enzyme_name || '') + '<br><b>Family:</b> ' + (e.enzyme_family || '') + '<br><b>Pathway:</b> ' + (e.pathway || '') + '<br><b>Substrate:</b> ' + (e.substrate_type || '') + '<br><br><b>Description:</b><br>' + (e.description || '') : '';
        const links = NETWORK.filter(n => {
          const enzyme = ENZYMES.find(x => x.enzyme_id === n.enzyme_id);
          return enzyme && enzyme.gene_symbol === sym;
        });
        const rows = links.map(link => {
          const g = GLYCANS.find(x => x.glycan_id === link.glycan_id);
          return '<tr class="border-b"><td class="px-4 py-2">' + (g ? g.snfg_name : '') + '</td><td class="px-4 py-2">' + (g ? g.composition : '') + '</td><td class="px-4 py-2">' + (link.linkage || '') + '</td><td class="px-4 py-2"><span class="px-2 py-1 bg-blue-100 text-blue-800 rounded text-xs">' + (link.evidence || '') + '</span></td></tr>';
        }).join('');
        document.getElementById('table-enzyme-glycans').innerHTML =
          '<table class="min-w-full text-left text-sm"><thead class="bg-gray-100"><tr><th class="px-4 py-2">Glycan</th><th class="px-4 py-2">Composition</th><th class="px-4 py-2">Linkage</th><th class="px-4 py-2">Evidence</th></tr></thead><tbody>' + rows + '</tbody></table>';
      }
      updateEnzymeView();
    </script>
  `;
  return pageTemplate('Enzyme Network', content, '/enzymes');
}

// ========================================================================
// TCGA Expression Page
// ========================================================================
function tcgaPage() {
  const tcga = DB.tcga_data;

  const exprChart = JSON.stringify({
    xAxis: { type: 'category', data: tcga.map(d => d.gene), axisLabel: { rotate: 45 } },
    yAxis: { type: 'value', name: 'Median RPKM' },
    series: [
      { name: 'Normal', type: 'bar', data: tcga.map(d => d.normal_rpkm), itemStyle: { color: '#00a65a' } },
      { name: 'PTC', type: 'bar', data: tcga.map(d => d.ptc_rpkm), itemStyle: { color: '#dd4b39' } }
    ],
    tooltip: { trigger: 'axis' },
    legend: { data: ['Normal', 'PTC'] }
  });

  const fcChart = JSON.stringify({
    xAxis: { type: 'category', data: tcga.map(d => d.gene), axisLabel: { rotate: 45 } },
    yAxis: { type: 'value', name: 'Fold Change (PTC/Normal)' },
    series: [{
      data: tcga.map(d => ({
        value: d.fold_change,
        itemStyle: { color: d.significant ? (d.fold_change > 1 ? '#dd4b39' : '#00a65a') : '#bbbbbb' }
      })),
      type: 'bar',
      label: { show: true, position: 'top', formatter: '{c}' }
    }],
    tooltip: { trigger: 'axis' }
  });

  const rows = tcga.map(g => `
    <tr class="border-b hover:bg-gray-50">
      <td class="px-4 py-2 text-sm font-medium">${g.gene}</td>
      <td class="px-4 py-2 text-sm">${g.family || ''}</td>
      <td class="px-4 py-2 text-sm">${g.normal_rpkm}</td>
      <td class="px-4 py-2 text-sm">${g.ptc_rpkm}</td>
      <td class="px-4 py-2 text-sm">${g.fold_change ? g.fold_change.toFixed(2) : 'N/A'}</td>
      <td class="px-4 py-2 text-sm">${g.p_value || ''}</td>
      <td class="px-4 py-2 text-sm">${g.significant ? '<span class="text-red-600 font-bold">*</span>' : 'ns'}</td>
    </tr>
  `).join('');

  const content = `
    <div class="mb-6">
      <h1 class="text-2xl font-bold text-gray-900">TCGA Glycogene Expression</h1>
      <p class="text-gray-600 mt-1">Real TCGA-THCA RNAseq data from Bones et al. 2018 (n=20 paired)</p>
    </div>
    <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
      <div class="bg-white rounded-lg shadow p-4">
        <h3 class="font-semibold text-gray-800 mb-2">Median RPKM by Group</h3>
        <div id="chart-tcga-expr" class="chart-container-lg"></div>
      </div>
      <div class="bg-white rounded-lg shadow p-4">
        <h3 class="font-semibold text-gray-800 mb-2">Fold Change & Significance</h3>
        <div id="chart-tcga-fc" class="chart-container-lg"></div>
      </div>
    </div>
    <div class="bg-white rounded-lg shadow overflow-x-auto">
      <table class="min-w-full text-left">
        <thead class="bg-gray-100 text-gray-700 text-xs uppercase">
          <tr><th class="px-4 py-2">Gene</th><th class="px-4 py-2">Family</th><th class="px-4 py-2">Normal RPKM</th><th class="px-4 py-2">PTC RPKM</th><th class="px-4 py-2">FC</th><th class="px-4 py-2">p-value</th><th class="px-4 py-2">Sig</th></tr>
        </thead>
        <tbody class="text-sm">${rows}</tbody>
      </table>
    </div>
    <script>
      echarts.init(document.getElementById('chart-tcga-expr')).setOption(${exprChart});
      echarts.init(document.getElementById('chart-tcga-fc')).setOption(${fcChart});
    </script>
  `;
  return pageTemplate('TCGA Expression', content, '/tcga');
}

// ========================================================================
// Diagnostic Tool Page
// ========================================================================
function diagnosticPage() {
  const lit = DB.lit_stats;
  const bn = lit.find(x => x.study === 'Zhang_2021' && x.variable && x.variable.includes('BN'));
  const ratio = lit.find(x => x.study === 'Kudelka_2023' && x.variable && x.variable.includes('G0F:G1F'));

  const bnRef = bn ? `HC: ${bn.group1_mean || 'N/A'} ± ${bn.group1_sd || 'N/A'}% | TC: ${bn.group2_mean || 'N/A'} ± ${bn.group2_sd || 'N/A'}% | AUC: ${bn.auc || 'N/A'}` : 'N/A';
  const ratioRef = ratio ? `HC median: ${ratio.group1_median || 'N/A'} | Recurrent: ${ratio.group2_median || 'N/A'} | AUC: ${ratio.auc || 'N/A'}` : 'N/A';

  const content = `
    <div class="mb-6">
      <h1 class="text-2xl font-bold text-gray-900">Diagnostic Tool</h1>
      <p class="text-gray-600 mt-1">Interactive calculators based on published biomarker performance</p>
    </div>
    <div class="grid grid-cols-1 md:grid-cols-2 gap-6">
      <div class="bg-white rounded-lg shadow p-6">
        <h3 class="font-semibold text-gray-800 mb-2">IgG BN Diagnostic Score (Zhang 2021)</h3>
        <p class="text-xs text-gray-500 mb-4">${bnRef}</p>
        <div class="space-y-3">
          <div><label class="text-sm">H3N5F1 (%)</label><input type="number" id="d-h3n5f1" value="15" class="w-full border rounded px-2 py-1 text-sm"></div>
          <div><label class="text-sm">H4N5F1 (%)</label><input type="number" id="d-h4n5f1" value="20" class="w-full border rounded px-2 py-1 text-sm"></div>
          <div><label class="text-sm">H5N5F1 (%)</label><input type="number" id="d-h5n5f1" value="10" class="w-full border rounded px-2 py-1 text-sm"></div>
          <div><label class="text-sm">H4N4F1 (%) Reference</label><input type="number" id="d-h4n4f1" value="55" class="w-full border rounded px-2 py-1 text-sm"></div>
          <button onclick="calcBN()" class="bg-blue-600 text-white px-4 py-2 rounded text-sm hover:bg-blue-700">Calculate BN Score</button>
        </div>
        <div id="d-result" class="mt-4 p-3 bg-gray-50 rounded text-sm font-mono"></div>
      </div>
      <div class="bg-white rounded-lg shadow p-6">
        <h3 class="font-semibold text-gray-800 mb-2">Serum Recurrence Predictor (Kudelka 2023)</h3>
        <p class="text-xs text-gray-500 mb-4">${ratioRef}</p>
        <div class="space-y-3">
          <div><label class="text-sm">G0F Intensity</label><input type="number" id="d-g0f" value="35" class="w-full border rounded px-2 py-1 text-sm"></div>
          <div><label class="text-sm">G1F Intensity</label><input type="number" id="d-g1f" value="50" class="w-full border rounded px-2 py-1 text-sm"></div>
          <button onclick="calcRatio()" class="bg-yellow-600 text-white px-4 py-2 rounded text-sm hover:bg-yellow-700">Calculate Ratio</button>
        </div>
        <div id="r-result" class="mt-4 p-3 bg-gray-50 rounded text-sm font-mono"></div>
      </div>
    </div>
    <script>
      function calcBN() {
        const h3 = parseFloat(document.getElementById('d-h3n5f1').value) || 0;
        const h4 = parseFloat(document.getElementById('d-h4n5f1').value) || 0;
        const h5 = parseFloat(document.getElementById('d-h5n5f1').value) || 0;
        const ref = parseFloat(document.getElementById('d-h4n4f1').value) || 0;
        const total = h3 + h4 + h5 + ref;
        const bn = total > 0 ? (h3 + h4 + h5) / total : 0;
        document.getElementById('d-result').innerText = 'BN Ratio: ' + bn.toFixed(4) + '\\n' + (bn > 0.35 ? 'BN ELEVATED (consistent with TC)' : 'BN not elevated (consistent with HC)');
      }
      function calcRatio() {
        const g0f = parseFloat(document.getElementById('d-g0f').value) || 0;
        const g1f = parseFloat(document.getElementById('d-g1f').value) || 0;
        if (g1f === 0) { document.getElementById('r-result').innerText = 'Error: G1F cannot be zero'; return; }
        const ratio = g0f / g1f;
        let msg = 'G0F:G1F = ' + ratio.toFixed(3) + '\\n';
        if (ratio > 0.73) msg += 'HIGH recurrence risk (ratio > 0.73)';
        else if (ratio > 0.53) msg += 'MODERATE recurrence risk (0.53-0.73)';
        else msg += 'LOW recurrence risk (ratio < 0.53)';
        document.getElementById('r-result').innerText = msg;
      }
    </script>
  `;
  return pageTemplate('Diagnostic Tool', content, '/diagnostic');
}

// ========================================================================
// Nomogram Page
// ========================================================================
function nomogramPage() {
  const lit = DB.lit_stats;
  const ca4 = lit.find(x => x.study === 'PTMC_2022_LNM' && x.variable === 'CA4');
  const a2 = lit.find(x => x.study === 'PTMC_2022_LNM' && x.variable && x.variable.includes('A2F0S0G'));

  const ca4Ref = ca4 ? `NLNM median: ${ca4.group1_median || 'N/A'} | LNM median: ${ca4.group2_median || 'N/A'}` : 'N/A';
  const a2Ref = a2 ? `NLNM median: ${a2.group1_median || 'N/A'} | LNM median: ${a2.group2_median || 'N/A'}` : 'N/A';

  const content = `
    <div class="mb-6">
      <h1 class="text-2xl font-bold text-gray-900">PTMC LNM Risk Nomogram</h1>
      <p class="text-gray-600 mt-1">Demonstration model based on PTMC nomogram (Front Oncol 2022). <b>Not for clinical use without validation.</b></p>
    </div>
    <div class="grid grid-cols-1 md:grid-cols-2 gap-6">
      <div class="bg-white rounded-lg shadow p-6">
        <h3 class="font-semibold text-gray-800 mb-2">Input Values</h3>
        <p class="text-xs text-gray-500 mb-4">CA4: ${ca4Ref}<br>A2F0S0G: ${a2Ref}</p>
        <div class="space-y-3">
          <div><label class="text-sm">CA4 (Tetraantennary complex)</label><input type="number" id="n-ca4" value="0.022" step="0.001" class="w-full border rounded px-2 py-1 text-sm"></div>
          <div><label class="text-sm">A2F0S0G (Galactosylation index)</label><input type="number" id="n-a2" value="0.60" step="0.001" class="w-full border rounded px-2 py-1 text-sm"></div>
          <button id="btn-calc-nom" class="bg-blue-600 text-white px-4 py-2 rounded text-sm hover:bg-blue-700">Calculate Risk</button>
        </div>
      </div>
      <div class="bg-white rounded-lg shadow p-6">
        <h3 class="font-semibold text-gray-800 mb-2">Risk Probability</h3>
        <div id="n-gauge" class="chart-container"></div>
        <div id="n-text" class="mt-2 text-center text-sm font-medium"></div>
      </div>
    </div>
    <script>
      let nomogramChart;
      function calcNom() {
        try {
          const ca4 = parseFloat(document.getElementById('n-ca4').value) || 0;
          const a2 = parseFloat(document.getElementById('n-a2').value) || 0;
          const logit = -2.5 + 3.2 * ca4 - 2.0 * a2;
          const prob = 1 / (1 + Math.exp(-logit));
          const pct = (prob * 100).toFixed(1);
          let level = prob < 0.3 ? 'Low' : prob < 0.6 ? 'Intermediate' : 'High';
          let color = prob < 0.3 ? '#00a65a' : prob < 0.6 ? '#f39c12' : '#dd4b39';
          document.getElementById('n-text').innerHTML = '<span style="color:' + color + '">' + level + ' Risk: ' + pct + '%</span>';
          if (typeof echarts === 'undefined') {
            document.getElementById('n-text').innerHTML += '<br><span style="color:#999;font-size:12px">Chart library not loaded</span>';
            return;
          }
          const gaugeEl = document.getElementById('n-gauge');
          if (!nomogramChart) {
            nomogramChart = echarts.init(gaugeEl);
          }
          nomogramChart.setOption({
            series: [{
              type: 'gauge',
              min: 0, max: 100,
              axisLine: { lineStyle: { width: 20, color: [[0.3, '#00a65a'], [0.6, '#f39c12'], [1, '#dd4b39']] } },
              pointer: { length: '60%', width: 5 },
              detail: { fontSize: 20, formatter: '{value}%' },
              data: [{ value: parseFloat(pct), name: 'Risk %' }]
            }]
          });
        } catch (e) {
          document.getElementById('n-text').innerHTML = '<span style="color:#dd4b39">Error: ' + e.message + '</span>';
          console.error(e);
        }
      }
      document.getElementById('btn-calc-nom').addEventListener('click', calcNom);
      calcNom();
    </script>
  `;
  return pageTemplate('Nomogram', content, '/nomogram');
}

// ========================================================================
// API Handlers
// ========================================================================
function apiResponse(data, status = 200) {
  return new Response(JSON.stringify(data), {
    status: status,
    headers: { 'Content-Type': 'application/json', 'Access-Control-Allow-Origin': '*' }
  });
}

function handleApi(path) {
  if (path === '/api/studies') return apiResponse(DB.studies);
  if (path === '/api/biomarkers') return apiResponse(DB.biomarkers);
  if (path === '/api/glycans') return apiResponse(DB.glycans);
  if (path === '/api/enzymes') return apiResponse(DB.enzymes);
  if (path === '/api/tcga') return apiResponse(DB.tcga_data);
  if (path === '/api/stats') return apiResponse(DB.stats);
  if (path === '/api/literature-stats') return apiResponse(DB.lit_stats);
  return apiResponse({ error: 'Not found' }, 404);
}

// ========================================================================
// Main Router
// ========================================================================
export default {
  async fetch(request, env, ctx) {
    const url = new URL(request.url);
    const path = url.pathname;

    if (path.startsWith('/api/')) {
      return handleApi(path);
    }

    switch (path) {
      case '/':
      case '/overview':
        return new Response(overviewPage(), { headers: { 'Content-Type': 'text/html; charset=utf-8' } });
      case '/literature':
        return new Response(literaturePage(), { headers: { 'Content-Type': 'text/html; charset=utf-8' } });
      case '/biomarkers':
        return new Response(biomarkersPage(), { headers: { 'Content-Type': 'text/html; charset=utf-8' } });
      case '/glycans':
        return new Response(glycansPage(), { headers: { 'Content-Type': 'text/html; charset=utf-8' } });
      case '/enzymes':
        return new Response(enzymesPage(), { headers: { 'Content-Type': 'text/html; charset=utf-8' } });
      case '/tcga':
        return new Response(tcgaPage(), { headers: { 'Content-Type': 'text/html; charset=utf-8' } });
      case '/diagnostic':
        return new Response(diagnosticPage(), { headers: { 'Content-Type': 'text/html; charset=utf-8' } });
      case '/nomogram':
        return new Response(nomogramPage(), { headers: { 'Content-Type': 'text/html; charset=utf-8' } });
      default:
        return new Response(overviewPage(), { headers: { 'Content-Type': 'text/html; charset=utf-8' } });
    }
  }
};
