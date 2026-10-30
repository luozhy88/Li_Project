# -*- coding: utf-8 -*-
"""Build the English structured "Materials and Methods / Results / Conclusion"
submission reference (word/character counts annotated per section)."""
from docx import Document
from docx.shared import Pt
from docx.oxml.ns import qn

METHODS = (
    "We used network pharmacology followed by validation in independent transcriptome cohorts. Ingredients of "
    "the core herbs selected from Xiao Ying Tang (XYT; 14 herbs including E Zhu, Bai Shao, Dang Gui, Chai Hu, "
    "Xia Ku Cao and Zhe Bei Mu) were taken from TCMSP and mapped to PubChem identifiers. Compound-associated "
    "genes came from BindingDB, BATMAN-TCM (prediction score >0.80) and the IUPHAR/BPS Guide to PHARMACOLOGY; "
    "intersecting them with 12,987 GeneCards thyroid-tumour-associated genes gave the XYT, Zhe Bei Mu (ZBM) and "
    "Xia Ku Cao (XKC) candidate sets, which were analysed with clusterProfiler for GO and KEGG enrichment. "
    "Expression was studied in two paired papillary thyroid carcinoma (PTC) cohorts, GSE60542 (24 patients) and "
    "GSE33630 (44 pairs), with paired limma models; differential expression required Benjamini\u2013Hochberg-"
    "adjusted P <0.05 and |log2 fold change| \u22651. Candidate-set enrichment was tested with one-sided Fisher "
    "tests against all measured genes and against the disease-associated subset, and cross-cohort agreement was "
    "checked; GSE33630 was also analysed separately by radiation exposure. \u03b2-Sitosterol, the compound shared "
    "by all three sets, was docked to human CYP3A4 (6UNI) and CaV3.2 (9AYG) as an exploratory analysis."
)

RESULTS = (
    "The disease intersections contained 78 ZBM, 842 XKC and 1,595 XYT genes. GO and KEGG enrichment centred on "
    "oxidative stress, calcium signalling, MAPK signalling and the thyroid cancer pathway. Each cohort measured "
    "20,848 genes; 666 in GSE60542 and 745 in GSE33630 met the differential-expression criteria. XYT and XKC "
    "were enriched in both cohorts, with odds ratios against disease-associated backgrounds of 1.80 and 2.08 in "
    "GSE60542 and 1.59 and 1.80 in GSE33630 (all adjusted P <0.001); ZBM was not enriched. Of the 666 genes "
    "found in GSE60542, 664 changed in the same direction in GSE33630 and 536 met the same criteria again, "
    "including 83 XYT candidates \u2014 among them the thyroid-differentiation markers DIO1 and TPO, both "
    "down-regulated, and LRRK2, DPP4 and FN1, all up-regulated. Enrichment held in both radiation-exposure "
    "strata. POLB, CACNA1H and CYP3A4 met the combined criteria in neither cohort. The best docking scores for "
    "\u03b2-sitosterol were \u22129.12 kcal/mol (CYP3A4) and \u22129.68 kcal/mol (CaV3.2), which we treat as "
    "exploratory."
)

DISCUSSION = (
    "This study ties the selected XYT core herbs to thyroid-disease genes and then finds those genes altered in "
    "paired PTC tissue. The enriched annotations \u2014 oxidative stress, calcium regulation, signalling \u2014 fit "
    "what is known of thyroid cancer biology, but they describe the candidate sets as a whole, not the activity "
    "of any single pathway. Because the three sets are nested, the wider XYT coverage is expected and cannot be "
    "read as synergy, and the non-significant ZBM result proves neither inactivity nor a difference from the "
    "other sets, resting as it does on only 77 measured genes. The replicated single-gene candidates were headed "
    "by genes already characteristic of PTC, including the down-regulated differentiation markers DIO1 and TPO. "
    "POLB fell slightly in both cohorts, but transcript abundance is not protein activity, and the docking poses "
    "add little: CYP3A4 lacked its haem cofactor, and the automatic search regions were never validated as "
    "binding sites. Limitations include target assignments that mix recorded and predicted associations, a broad "
    "disease list that favours well-studied genes, and the young GSE33630 patients. Cross-cohort agreement "
    "supports reproducible disease associations, not a response to herbal treatment. Follow-up work should bring "
    "together concentration\u2013response experiments, non-malignant controls, target perturbation and binding or "
    "functional assays."
)

CONCLUSION = (
    "XYT and XKC candidate genes were over-represented among genes altered in paired PTC tissue in two "
    "independent cohorts, including against disease-associated backgrounds. Eighty-three XYT candidates, among "
    "them the thyroid-differentiation markers DIO1 and TPO, met the differential-expression criteria in both "
    "cohorts and are the genes to take forward into functional work. ZBM showed no clear enrichment, and none "
    "of the three preselected genes (POLB, CACNA1H, CYP3A4) met the combined criteria. These are disease "
    "associations, not evidence of treatment effect; whether the selected herbs have any therapeutic role in "
    "thyroid cancer still needs testing in treated models with chemically characterised preparations."
)

doc = Document()
style = doc.styles["Normal"]
style.font.name = "Times New Roman"
style.font.size = Pt(11)
style.element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")

doc.add_heading(
    "Network pharmacology and independent transcriptome validation of core herbs selected from "
    "Xiao Ying Tang in thyroid cancer", level=1)
doc.add_paragraph(
    'Structured "Materials and Methods / Results / Discussion / Conclusion" for submission (English; word '
    "and character counts annotated per section).")

for title, text in [("Materials and Methods", METHODS), ("Results", RESULTS),
                    ("Discussion", DISCUSSION), ("Conclusion", CONCLUSION)]:
    doc.add_heading(title, level=2)
    doc.add_paragraph(text)
    doc.add_paragraph(f"({len(text.split())} words, {len(text)} characters)")

out = "output/20261027_27_before_final_review/08_Structured_sections_EN_reference.docx"
doc.save(out)
print("saved:", out)
for title, text in [("Materials and Methods", METHODS), ("Results", RESULTS),
                    ("Discussion", DISCUSSION), ("Conclusion", CONCLUSION)]:
    print(f"{title}: {len(text.split())} words, {len(text)} chars")
