# Response to Reviewers

## Manuscript: 针灸穴位-基因关联数据库的构建及数据挖掘分析

---

We would like to express our sincere gratitude to the Editor and Reviewers for their thorough and constructive comments on our manuscript. We have carefully revised the manuscript according to the suggestions, and have also performed a comprehensive data quality audit and cleanup on the underlying database. Below, we provide point-by-point responses to each major concern raised.

---

## 1. Major Concern: Gene Entity Recognition Errors

**Reviewer Comment:** The database contains numerous non-gene entities misidentified as gene symbols (e.g., METABOLOMIC, POSTOPERATIVE, STROKE, SYNDROME, WEB, T2DM, BBB, etc.). In addition, gene family names (AKT, PI3K, MAPK) were treated as single genes, and unofficial symbols (LC3, BECLIN) were not mapped to official HGNC symbols.

**Response:** We fully acknowledge this critical data quality issue and have conducted a comprehensive cleanup:

1. **Removed non-gene entities:** We established a dictionary-based filtering mechanism and removed 73 non-gene entities from the gene table, including technical terms (METABOLOMIC, CYTOSCAPE), disease names (STROKE, SYNDROME, T2DM), clinical abbreviations (WOMAC, HOMA, NIHSS), and ordinary English words (RESPONSE, INSIGHTS, OVERLAPPING). All associations involving these entities have been deleted from `article_gene`.

2. **Annotated gene families/pathways:** The 10 gene families/pathways (AKT, PI3K, MAPK, WNT, JAK, ERK, MMP, BCL, P38, STAT) are now annotated separately in the database with `entity_type = 'gene_family'` and are excluded from single-gene statistics. We have added a dedicated section (Section 3.3) to report these families separately.

3. **Mapped unofficial symbols:** LC3 has been mapped to MAP1LC3A/B/C, and BECLIN/BECLIN1 have been mapped to BECN1, following HGNC nomenclature.

4. **Revised statistics:** The corrected database now contains **1,997 valid single genes** (previously claimed 2,083), plus **10 gene families/pathways**.

**Revisions in manuscript:**
- Section 2.3: Revised gene standardization workflow to include dictionary filtering, family annotation, and HGNC mapping.
- Section 3.3: Split into "High-frequency single genes" and "Gene families/pathways".
- Section 4 (Limitations): Added transparent disclosure of these data quality issues.

---

## 2. Major Concern: Co-occurrence Overinterpreted as "Regulation"

**Reviewer Comment:** Co-occurrence in abstracts does not establish regulatory causation. The manuscript repeatedly interprets co-occurrence frequency as "association strength" and even "regulatory specificity," which constitutes methodological overinterpretation. Additionally, acupoint combination effects cannot be attributed to single acupoints.

**Response:** We agree completely and have softened all causal language throughout the manuscript:

1. **Abstract & Conclusion:** Changed "acupoint regulation of gene expression" to "acupoint-gene co-occurrence patterns" and "research attention."

2. **Section 3.5:** Added explicit methodological disclaimer: "Co-occurrence frequency reflects research attention, not direct regulatory causation. Genes may be mentioned as background mechanisms, disease biomarkers, or controls."

3. **Acupoint combination caveat:** Added explanation that many studies use acupoint combinations (e.g., Baihui + Taichong, Zusanli + Sanyinjiao), and co-occurrence analysis cannot disentangle combination effects into single-acupoint contributions. The conclusion "relative specificity" has been revised to "differences in research attention across acupoints."

4. **Discussion:** Added statement that "relative specificity" conclusions require controlled experiments (acupoint vs. non-acupoint, meridian vs. non-meridian) for validation.

**Revisions in manuscript:**
- Abstract, Sections 3.5, 4, and 5: All causal/regulatory language softened to co-occurrence/attention language.

---

## 3. Major Concern: Year Data Anomaly (2025–2026)

**Reviewer Comment:** The 17-fold jump from 19 articles in 2024 to 111 in the first 4 months of 2026 likely reflects PubMed "ahead of print" artifacts rather than genuine research growth. No technical explanation was provided.

**Response:** We have added a detailed technical explanation:

1. **Section 3.1:** Added note: "The 2025–2026 data contain a substantial number of Online ahead of print articles, whose official publication years may be delayed. Therefore, the surge in 2025–2026 may partly reflect the proliferation of online-first publication models, rather than a pure increase in literature volume."

2. **Database:** Added `pub_status` field to the `articles` table to distinguish "formally published" from "ahead of print."

3. **Trend analysis:** We caution readers against over-interpreting the 2025–2026 trend until final publication years are confirmed.

**Revisions in manuscript:**
- Section 3.1: Added technical explanation and caveat.
- Section 4 (Limitations): Added discussion of year-label instability.

---

## 4. Important Concern: Inclusion/Exclusion Criteria Contradictions

**Reviewer Comment:** (a) 74 Chinese-language articles contradict the "English-only" inclusion criterion; (b) 36 Review/Letter articles contradict the "exclude reviews" criterion.

**Response:** We have revised the criteria and marked problematic records:

1. **Language issue:** We clarified that the 74 "chi" records are English articles published in Chinese journals (e.g., *Zhen ci yan jiu*, *Zhongguo zhen jiu*) with English online editions. The inclusion criterion has been revised to "English (including English online editions of Chinese journals)."

2. **Article type issue:** We acknowledge that 36 Review/Letter articles were inadvertently retained. These have been flagged in the database with `data_quality_note` and will be excluded in the next data release. The exclusion criteria have been expanded to explicitly exclude: Review, Systematic Review, Letter, Commentary, and Editorial.

3. **Expanded exclusion criteria:** Added items (4)–(6): animal studies without species/model reporting; cell studies without stimulation method; ahead-of-print without complete abstracts.

**Revisions in manuscript:**
- Section 2.2: Revised inclusion/exclusion criteria with expanded detail.
- Section 4 (Limitations): Disclosed the language and article-type issues.

---

## 5. Important Concern: Disease Classification Confusion

**Reviewer Comment:** Inflammation, pain, and stress are pathophysiological processes or symptoms, not diseases. Combining them with true diseases (tumor, stroke) creates taxonomic confusion. Inflammatory/inflammation were counted separately.

**Response:** We have implemented a two-tier classification system:

1. **Merged synonyms:** `inflammatory` and `inflammation` have been merged into a single entity (`inflammation`, 157 articles after merge).

2. **Two-tier taxonomy:**
   - **Pathological process/symptom:** inflammation (157), pain (69), stress (59), withdrawal (25), depression (20), anxiety (13), etc.
   - **Disease:** tumor (38), stroke (31), arthritis (30), cancer (17), osteoarthritis (13), rheumatoid arthritis (12), etc.

3. **Section 3.4:** Revised to report diseases and pathological processes in separate subsections.

**Revisions in manuscript:**
- Section 2.3: Added two-tier classification strategy.
- Section 3.4: Restructured into disease-class and pathological-process-class statistics.

---

## 6. Suggestion: Insufficient Analytical Depth

**Reviewer Comment:** As a "database + data mining" paper, only frequency counting and co-occurrence are insufficient. Network topology, enrichment analysis, and statistical testing are expected.

**Response:** We have addressed this in two ways:

1. **Current analysis:** We acknowledge the current analytical depth is limited and have added explicit suggestions in Section 2.4 for:
   - Network topology analysis (degree centrality, betweenness centrality)
   - KEGG/GO enrichment analysis (clusterProfiler/DAVID)
   - Hypergeometric test / chi-square test for co-occurrence significance

2. **Future work:** These deeper analyses are planned for the next database release and will be reported in a follow-up study. We believe the current manuscript provides a solid foundation (database construction + basic descriptive statistics) and the proposed methods offer a clear roadmap.

**Revisions in manuscript:**
- Section 2.4: Added suggested analytical methods.
- Section 4 (Limitations): Acknowledged analytical depth limitation.
- Section 4 (Future work): Committed to implementing these methods.

---

## 7. Suggestion: "First" Claim Risk

**Reviewer Comment:** The claim "首个" (first) is high-risk. If reviewers find existing similar databases, the paper may be rejected outright.

**Response:** We have completely removed the "首个" claim. The revised introduction states: "Currently, although some acupuncture-related databases exist, a literature integration platform specifically focused on acupoint-gene associations using text mining methods is still lacking." This formulation highlights our methodological contribution (text mining + systematic integration) rather than claiming absolute priority.

**Revisions in manuscript:**
- Section 1 (Introduction): Removed "首个" and "填补空白" language.
- Section 4 (Discussion): Revised to emphasize methodological features (text mining, GB/T 12346—2021 standardization, 2010–2026 coverage).

---

## 8. Additional Revisions

### 8.1 Acupoint Synonym Merging
We confirmed that the database-level `acupoints` table correctly merges synonyms (e.g., ST36/Zusanli/Sanli → 足三里, 109 articles; Baihui/GV20 → 百会, 26 articles). We have verified that all statistics in the revised manuscript are based on the merged table, not the outdated `summary.json`.

### 8.2 Data Quality Transparency
To demonstrate our commitment to data quality, we have:
- Added a `data_quality_note` field to 261 articles flagging language, type, or year issues.
- Generated `problematic_genes.csv` listing all 86 problematic entities with correction suggestions.
- Generated `flagged_articles.csv` listing all 261 flagged articles with PMIDs and quality notes.
- These files are available as supplementary materials upon request.

### 8.3 Figures Updated
All figures have been regenerated based on the cleaned data:
- **Fig 1:** Dashboard now shows "Valid Genes: 1,997" (previously 2,083).
- **Fig 2:** Top genes now exclude AKT, PI3K, MAPK families; BECN1 and MAP1LC3A/B/C appear with official symbols.
- **Fig 4:** Network graph now excludes non-gene entities and gene families.

---

## Summary of Changes

| Issue | Action Taken | Manuscript Section |
|-------|-------------|-------------------|
| Non-gene entities (73) | Removed from gene statistics | 2.3, 3.3, 4 |
| Gene families (10) | Annotated separately | 2.3, 3.3 |
| Unofficial symbols (3) | Mapped to HGNC | 2.3, 3.3 |
| Co-occurrence = regulation | Softened to "attention" throughout | Abstract, 3.5, 4, 5 |
| Year anomaly | Added technical explanation | 3.1, 4 |
| Language contradiction | Clarified criteria | 2.2, 4 |
| Article type contradiction | Expanded exclusion criteria | 2.2, 4 |
| Disease classification | Two-tier taxonomy + merge synonyms | 2.3, 3.4 |
| "First" claim | Removed, emphasized methodology | 1, 4 |
| Analytical depth | Added suggested methods | 2.4, 4 |
| Data transparency | Generated quality reports | Supplementary |

---

We believe the revised manuscript now addresses all major concerns regarding data quality, methodological rigor, and interpretive caution. We hope the Editor and Reviewers will find the revisions satisfactory and consider the manuscript for publication.

**Corresponding Author**
[Name]
[Affiliation]
[Email]
[Date]
