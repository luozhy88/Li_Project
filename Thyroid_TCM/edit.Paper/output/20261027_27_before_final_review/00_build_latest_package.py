from pathlib import Path
from zipfile import ZipFile
import shutil, subprocess, tempfile, hashlib
import pandas as pd
from openpyxl import load_workbook
from openpyxl.styles import Font, PatternFill
from openpyxl.utils import get_column_letter

root = Path('/home/data1/zhiyu/data/Github/Li_Project/Thyroid_TCM/edit.Paper/output')
out = root / '20261027_26_latest_submission'
out.mkdir(parents=True, exist_ok=True)
(out/'09_figures').mkdir(exist_ok=True)
(out/'10_editable_sources').mkdir(exist_ok=True)

src25 = root/'20261025_25_independent_cohort'
src19 = root/'20261025_19_JCPSP_submission'

# The package's own 10_editable_sources manuscripts are now authoritative
# (workflow figure added, figures renumbered, package-local image paths).
for name in ['01_Manuscript_EN.md', '02_Manuscript_CN.md']:
    shutil.copy2(out/'10_editable_sources'/name, out/name)

pandoc = '/home/zhiyu/miniconda3/envs/test/bin/pandoc'
for md, docx in [(out/'01_Manuscript_EN.md', out/'01_Manuscript_EN.docx'),
                 (out/'02_Manuscript_CN.md', out/'02_Manuscript_CN.docx'),
                 (out/'04_CCHTS_submission_guide.md', out/'04_CCHTS_submission_guide.docx')]:
    subprocess.run([pandoc, md.name, '--from', 'markdown+superscript', '-o', docx.name,
                    '--standalone'], check=True, cwd=out)

# Cover letters are maintained in the package (CCHTS version, EN + CN).
for cl in ['03_Cover_letter.md', '03_Cover_letter_CN.md']:
    shutil.copy2(out/'10_editable_sources'/cl, out/cl)
    subprocess.run([pandoc, cl, '-o', cl.replace('.md', '.docx'), '--standalone'], check=True, cwd=out)

# Figure 1 (technical route) and the graphical abstract are rendered from graphviz sources.
gv = out/'10_editable_sources/Figure_1_workflow.gv'
subprocess.run(['dot', '-Tpng', '-Gdpi=300', str(gv), '-o', str(out/'09_figures/Figure_1_workflow.png')], check=True)
subprocess.run(['dot', '-Tpdf', str(gv), '-o', str(out/'10_editable_sources/Figure_1_workflow.pdf')], check=True)
ga = out/'10_editable_sources/Graphical_Abstract.gv'
subprocess.run(['dot', '-Tpng', '-Gdpi=300', str(ga), '-o', str(out/'09_figures/Graphical_Abstract.png')], check=True)

# Main figures renumbered after inserting the workflow as Figure 1.
figs = [
    ('Figure_2_GO.png', src19/'09_figures/01_Figure_1.png'),
    ('Figure_3_KEGG.png', src19/'09_figures/02_Figure_2.png'),
    ('Figure_4_GSE60542.png', src19/'09_figures/03_Figure_3.png'),
    ('Supplementary_Figure_S1_gene_overlap.png', src19/'09_figures/04_Supplementary_Figure_S1.png'),
    ('Supplementary_Figure_S2_docking.png', src19/'09_figures/05_Supplementary_Figure_S2.png'),
    ('Supplementary_Figure_S3_GSE33630.png', src25/'16_Supplementary_Figure_S3.png'),
]
for name, path in figs: shutil.copy2(path, out/'09_figures'/name)

# Supplementary workbook: original supplementary ZIP is canonical; GSE33630
# analysis-directory CSVs are added only when their content is not already
# present (hash-based dedup; GSE33630_04 is additionally dropped because its
# values are a column subset of 28_GSE33630_set_enrichment).
tmp = Path(tempfile.mkdtemp(prefix='suppl_tables_'))
with ZipFile(root/'20261025_19_JCPSP_submission.zip') as z: z.extractall(tmp/'original')
nested = next((tmp/'original').rglob('05_Supplementary_tables.zip'))
with ZipFile(nested) as z: z.extractall(tmp/'tables')

entries = [
 ('01_target_counts.csv','T01_target_counts','Candidate-set sizes and disease-list overlap counts for ZBM, XKC and XYT','Results, first paragraph'),
 ('02_herb_inventory.csv','T02_herb_inventory','Herb inventory: unique PubChem CIDs per herb within each set','Materials and Methods; Results, first paragraph'),
 ('03_enrichment_results.csv','T03_enrichment_results','Full GO and KEGG enrichment results for the three candidate sets','Results; Figures 2-3'),
 ('04_figure2_source.csv','T04_figure2_source','Source data for Figure 2 (selected GO terms)','Figure 2'),
 ('05_figure3_source.csv','T05_figure3_source','Source data for Figure 3 (selected KEGG pathways)','Figure 3'),
 ('06_docking_audit.csv','T06_docking_audit','Docking audit: receptors, engines, best scores and search-box references','Materials and Methods; Supplementary Figure S2'),
 ('07_disease_overlap_genes_ZBM.csv','T07_overlap_genes_ZBM','ZBM disease-overlapping genes (n = 78)','Results, first paragraph'),
 ('08_disease_overlap_genes_XKC.csv','T08_overlap_genes_XKC','XKC disease-overlapping genes (n = 842)','Results, first paragraph'),
 ('09_disease_overlap_genes_XYT.csv','T09_overlap_genes_XYT','XYT disease-overlapping genes (n = 1,595)','Results, first paragraph'),
 ('10_gene_set_validation.csv','T10_gene_set_validation','Gene-set reconstruction validation checks','Materials and Methods'),
 ('11_target_source_counts.csv','T11_target_source_counts','Retained compound-gene records by target resource','Materials and Methods'),
 ('12_sample_manifest.csv','T12_GSE60542_sample_manifest','GSE60542 sample manifest with pairing and inclusion flags','Materials and Methods'),
 ('13_paired_samples.csv','T13_GSE60542_paired_samples','GSE60542 paired patient-tissue observations (n = 48)','Materials and Methods'),
 ('14_all_gene_differential_expression.csv','T14_GSE60542_all_gene_DE','GSE60542 all-gene paired limma results (20,848 genes)','Figure 4A'),
 ('15_candidate_gene_differential_expression.csv','T15_GSE60542_candidate_DE','GSE60542 differential-expression results for candidate-set genes','Results'),
 ('16_candidate_set_enrichment.csv','T16_GSE60542_set_enrichment','GSE60542 candidate-set enrichment (six one-sided Fisher tests)','Figure 4B'),
 ('17_candidate_coverage.csv','T17_GSE60542_coverage','GSE60542 candidate-set coverage on the measured background','Results'),
 ('18_probe_mapping.csv','T18_GSE60542_probe_mapping','GSE60542 probe-to-gene mapping summary','Materials and Methods'),
 ('19_gene_probe_counts.csv','T19_GSE60542_probe_counts','GSE60542 retained probe counts per gene','Materials and Methods'),
 ('23_transcriptome_enrichment_with_CI.csv','T20_GSE60542_enrichment_CI','GSE60542 candidate-set enrichment with exact 95% confidence intervals','Figure 4B'),
 ('25_GSE33630_sample_manifest.csv','T21_GSE33630_sample_manifest','GSE33630 sample manifest with pairing and inclusion flags','Materials and Methods'),
 ('26_GSE33630_verified_pairs.csv','T22_GSE33630_verified_pairs','GSE33630 verified matched pairs (n = 44)','Materials and Methods'),
 ('27_GSE33630_all_gene_DE.csv','T23_GSE33630_all_gene_DE','GSE33630 all-gene paired limma results (20,848 genes)','Supplementary Figure S3A'),
 ('28_GSE33630_set_enrichment.csv','T24_GSE33630_set_enrichment','GSE33630 candidate-set enrichment with exact 95% confidence intervals','Supplementary Figure S3B'),
 ('29_cross_cohort_gene_comparison.csv','T25_cross_cohort_genes','Cross-cohort per-gene comparison of limma results (20,848 genes)','Supplementary Figure S3C'),
 ('30_cross_cohort_summary.csv','T26_cross_cohort_summary','Cross-cohort replication summary (direction and repeated criteria)','Results'),
 ('31_exposure_sensitivity.csv','T27_exposure_sensitivity','GSE33630 exposure-stratified enrichment sensitivity analysis (12 tests)','Results'),
 ('32_exposed_stratum_DE.csv','T28_exposed_stratum_DE','GSE33630 exposed-stratum all-gene differential-expression results','Results'),
 ('33_unexposed_stratum_DE.csv','T29_unexposed_stratum_DE','GSE33630 unexposed-stratum all-gene differential-expression results','Results'),
 ('36_GSE33630_probe_QC.csv','T30_GSE33630_probe_QC','GSE33630 probe-to-gene mapping quality control','Materials and Methods'),
 ('37_GSE33630_model_QC.csv','T31_GSE33630_model_QC','GSE33630 paired-model quality control','Materials and Methods'),
 ('GSE33630_05_candidate_coverage.csv','T32_GSE33630_coverage','GSE33630 candidate-set coverage on the measured background','Results'),
 ('GSE33630_06_candidate_gene_DE.csv','T33_GSE33630_candidate_DE','GSE33630 differential-expression results for candidate-set genes','Results'),
 ('GSE33630_13_patient_group_samples.csv','T34_GSE33630_patient_groups','GSE33630 paired patient-tissue observations used in the model','Materials and Methods'),
 ('GSE33630_14_patient_group_expression.csv','T35_GSE33630_expression','GSE33630 averaged log2 expression matrix (88 observations x 20,848 genes)','Materials and Methods'),
 ('GSE33630_15_probe_mean_manual_check.csv','T36_GSE33630_probe_mean_QC','GSE33630 manual check of probe averaging','Materials and Methods'),
 ('GSE33630_17_design_matrix.csv','T37_GSE33630_design_matrix','GSE33630 paired-model design matrix','Materials and Methods'),
 ('GSE33630_19_expression_distribution.csv','T38_GSE33630_distribution','GSE33630 per-sample expression distribution summary','Materials and Methods'),
]

# T39: replicated XYT candidate genes, derived from the cross-cohort table
# and the candidate gene lists.
cc = pd.read_csv(tmp/'tables/29_cross_cohort_gene_comparison.csv')
xyt_genes = set(pd.read_csv(tmp/'tables/09_disease_overlap_genes_XYT.csv')['Com_Gene'])
xkc_genes = set(pd.read_csv(tmp/'tables/08_disease_overlap_genes_XKC.csv')['Com_Gene'])
zbm_genes = set(pd.read_csv(tmp/'tables/07_disease_overlap_genes_ZBM.csv')['Com_Gene'])
rep = cc[cc.DE_threshold_first & cc.DE_threshold_second & cc.same_direction].copy()
rep = rep[rep.gene.isin(xyt_genes)].copy()
rep['mean_abs'] = (rep.logFC_first.abs() + rep.logFC_second.abs())/2
rep = rep.sort_values('mean_abs', ascending=False)
t39 = pd.DataFrame({
    'gene': rep.gene,
    'in_XYT': True,
    'in_XKC': rep.gene.isin(xkc_genes),
    'in_ZBM': rep.gene.isin(zbm_genes),
    'log2FC_GSE60542': rep.logFC_first.round(6),
    'adj_P_GSE60542': rep['adj.P.Val_first'].round(6),
    'log2FC_GSE33630': rep.logFC_second.round(6),
    'adj_P_GSE33630': rep['adj.P.Val_second'].round(6),
})

xlsx = out/'05_Supplementary_Tables.xlsx'
with pd.ExcelWriter(xlsx, engine='openpyxl') as writer:
    readme_rows = []
    dfs = []
    for fname, sheet, desc, loc in entries:
        p = (tmp/'tables'/fname) if (tmp/'tables'/fname).exists() else src25/fname
        try: df = pd.read_csv(p, low_memory=False)
        except Exception: df = pd.read_csv(p, low_memory=False, encoding='latin1')
        dfs.append((sheet, df))
        readme_rows.append([sheet, desc, loc, fname, len(df)])
    dfs.append(('T39_replicated_candidates', t39))
    readme_rows.append(['T39_replicated_candidates',
                        'XYT candidate genes meeting the differential-expression criteria in both cohorts (n = 83)',
                        'Results',
                        'derived from 29_cross_cohort_gene_comparison.csv and candidate gene lists',
                        len(t39)])
    pd.DataFrame(readme_rows, columns=['Sheet','Description','Related manuscript content','Source CSV','Rows']).to_excel(writer, sheet_name='README', index=False)
    for sheet, df in dfs:
        df.to_excel(writer, sheet_name=sheet, index=False)
wb = load_workbook(xlsx)
for ws in wb.worksheets:
    ws.freeze_panes = 'A2'; ws.auto_filter.ref = ws.dimensions
    for c in ws[1]:
        c.font = Font(bold=True, color='FFFFFF'); c.fill = PatternFill('solid', fgColor='1F4E78')
    for col in ws.iter_cols():
        n = max([len(str(c.value)) if c.value is not None else 0 for c in col[:100]] or [10])
        ws.column_dimensions[get_column_letter(col[0].column)].width = min(42, max(10, n+2))
wb.save(xlsx)

for src_name, dst_name in [('21_manuscript_checks.json','21_manuscript_checks.json'),
                           ('22_独立队列验证结果与审核.md','22_GSE33630_validation_review.md'),
                           ('26_chinese_review_checks.json','26_chinese_review_checks.json'),
                           ('27_cover_letter_current.md','27_cover_letter_current.md')]:
    if (src25/src_name).exists(): shutil.copy2(src25/src_name, out/'10_editable_sources'/dst_name)
print('CREATED', out)
print('SHEETS', len(load_workbook(xlsx, read_only=True).sheetnames))
