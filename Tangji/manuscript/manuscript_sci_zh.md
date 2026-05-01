# ThyGlycoPortal: 基于真实文献数据的甲状腺癌糖基化交互分析平台构建与验证

> **建议英文题目**: ThyGlycoPortal: a literature-derived interactive platform for thyroid cancer glycomics integrating real published summary statistics

---

## 摘要 (Abstract)

**背景**：甲状腺癌（Thyroid Cancer, TC）是全球发病率增长最快的恶性肿瘤之一。蛋白质糖基化修饰在肿瘤发生、侵袭和转移中发挥关键调控作用，已成为液体活检新型标志物研究的前沿领域。然而，当前甲状腺癌糖组学研究成果呈高度碎片化分布，缺乏专门的知识整合平台，临床医生和研究人员难以系统性地获取、比较和转化这些发现。

**方法**：本研究系统检索PubMed/Medline数据库（截至2025年），按照系统评价与元分析优先报告规范（PRISMA）流程筛选纳入11项研究（其中9项提供原始糖组学数据的原创研究，2篇用于背景知识补充的综述）。建立双人独立提取+交叉验证的数据策展流程，从文献图表和表格中提取真实汇总统计量（均值±标准差、中位数、AUC、p值等），并建立严格的数据溯源与质量保证体系。采用SQLite关系型数据库（第三范式）整合糖链结构、临床分组、生物标志物、糖基转移酶网络和TCGA基因表达数据，基于R Shiny构建交互式Web平台。依据FAIR原则进行数据可发现性、可访问性、互操作性和可重用性自评估。

**结果**：ThyGlycoPortal v1.0涵盖9项提供原始数据的原创研究（2018–2025）、2篇背景综述、24种核心N-糖链结构、17个临床分组、7个验证/候选生物标志物、8个关键糖基转移酶和12条酶-糖链调控关联。平台收录了20对TCGA-THCA配对样本的糖基因真实median RPKM数据（Bones et al. 2018）。文献计量学分析显示，该领域在2018年后进入活跃期，MALDI-TOF MS占主导方法学（72.7%），样本类型以血清（45.5%）和组织（36.4%）为主。平台提供9大交互模块，支持3类典型应用场景：临床决策支持、跨研究糖链比较、多组学整合分析。FAIR自评总分为4.0/5.0。

**结论**：ThyGlycoPortal是首个专注于甲状腺癌的糖基化交互分析平台，填补了该领域专用知识库和临床转化工具的空白。平台所有数据均来源于已发表文献的真实汇总统计量或公开数据库（TCGA），为甲状腺癌糖组学知识的系统整合与临床转化提供了可靠的数据基础设施。

**关键词**：甲状腺癌；糖基化；糖组学；生物标志物；交互平台；数据库；真实数据；FAIR原则

---

## 1. 引言 (Introduction)

### 1.1 甲状腺癌的临床挑战与精准医学需求

甲状腺癌是内分泌系统最常见的恶性肿瘤，全球年龄标准化发病率在过去三十年中增长了约300%[1]。2020年全球新发病例约58.6万例，位列女性常见癌症第5位[2]。其中，乳头状甲状腺癌（Papillary Thyroid Carcinoma, PTC）约占全部甲状腺癌的85%以上，尽管总体10年生存率超过90%，但仍有10–20%的患者出现局部复发、颈部淋巴结转移（Lymph Node Metastasis, LNM）或远处转移，显著影响预后和生活质量[3,4]。

当前甲状腺癌的术前诊断主要依赖超声影像学检查和超声引导下细针穿刺活检（FNA）。然而，FNA细胞学诊断存在固有的局限性：约15–30%的结节被归类为"不确定"（Bethesda III–IV类），这部分患者往往需要通过诊断性手术切除来明确病理诊断，导致不必要的手术创伤和医疗负担[5,6]。因此，开发客观、无创、可重复的辅助诊断标志物，以及能够预测转移和复发风险的预后分层工具，是甲状腺癌精准医学领域的迫切需求[7]。

### 1.2 糖基化修饰在癌症生物学中的核心地位

蛋白质糖基化是真核生物中最普遍、最复杂的翻译后修饰之一，据估计超过50%的人类蛋白质带有糖链修饰[8]。糖基化广泛参与蛋白质折叠、细胞间识别、信号转导、免疫应答调控和肿瘤微环境塑造等关键生物学过程[9]。在肿瘤发生发展过程中，糖基化模式发生系统性重编程，表现为：

- **岩藻糖基化升高**：核心岩藻糖基化和外周岩藻糖基化在多种癌症中上调，与肿瘤侵袭性相关[10]；
- **半乳糖基化降低**：末端半乳糖残基减少（如G0F增加、G1F/G2F减少）是多种恶性肿瘤的共同特征[11]；
- **唾液酸化异常**：α2,3-和α2,6-连接唾液酸的比例改变影响肿瘤细胞的免疫逃逸和转移潜能[12]；
- **分支结构增多**：β1,6-GlcNAc分支（由MGAT5催化）与肿瘤的恶性程度和转移能力密切相关[13]。

近年来，基于血清或血浆N-糖组学的液体活检策略显示出巨大的临床转化潜力。通过高通量质谱技术（如MALDI-TOF MS）可以在微量体液样本中同时检测数十种糖链结构的相对丰度，构建高维度的"糖指纹"用于疾病诊断和预后评估[14,15]。

### 1.3 甲状腺癌糖组学研究的碎片化困境

尽管甲状腺癌糖组学研究在过去十年取得了显著进展（含9项原创研究及2篇领域综述），但核心实验数据呈高度碎片化分布。具体表现为：

**第一**，研究规模小而分散。已发表的原创研究多为单中心、小样本设计（样本量通常介于20–100例之间），不同研究采用的糖链命名系统、统计方法和报告标准各异，导致跨研究比较和Meta分析困难重重。

**第二**，数据可及性受限。原始糖组学数据（如质谱原始文件、个体水平糖链丰度矩阵）很少随论文公开发布，研究者通常仅报告汇总统计量（均值、标准差、AUC）和代表性图表。这限制了独立验证和二次分析的可能性。

**第三**，知识整合不足。糖链变化趋势、标志物性能参数、糖基转移酶表达数据和相关临床变量散落在不同文献中，临床医生和基础研究人员缺乏统一入口来系统性地浏览、查询和整合这些发现。

### 1.4 现有数据库与生物信息学工具的局限性

当前生物医学数据生态系统中存在多个与糖组学相关的数据库和分析工具，但它们在甲状腺癌领域的适用性存在显著缺陷（表1）：

**通用型糖组学数据库**如GlycoPOST[16]、UniCarbKB[17]和GlyConnect[18]主要聚焦于糖链结构注释、质谱原始数据存储或糖蛋白组学实验数据。它们虽然数据量庞大，但缺乏对特定癌种的深度整合，不提供临床标志物性能参数，也不支持交互式的跨研究比较分析。

**癌症多组学数据库**如TCGA[19]和THPA[20]虽覆盖甲状腺癌，但侧重点分别为基因组/转录组和蛋白质免疫组化图像。TCGA中糖基转移酶的表达信息深埋于庞大的表达矩阵中，缺乏针对性的注释和可视化工具；THPA则完全不涉及糖链结构层面的信息。

在数据来源方面，现有平台大多依赖作者自行上传的实验数据或从公共仓库批量下载的原始文件，对文献中已发表的汇总统计量（如均值、标准差、AUC）的系统性整合相对薄弱。对于临床医生和生物标志物研究者而言，汇总统计量往往比原始质谱数据更具直接参考价值，因为这些统计量已经过质量控制、归一化和统计检验，可直接用于标志物性能评估和临床决策建模。

### 1.5 研究目标

基于上述领域空白和技术挑战，本研究旨在构建ThyGlycoPortal——一个专注于甲状腺癌糖基化研究的文献驱动型交互分析平台。本研究的核心贡献包括：

1. **领域专精的知识整合**：系统收集并标准化2018–2025年间甲状腺癌糖组学文献的真实汇总统计量，建立首个该领域的专用关系型数据库；
2. **数据质量保证与溯源**：建立严格的数据策展流程和文献溯源机制，确保所有展示数据均可追溯至原始同行评议文献；
3. **交互式临床转化工具**：开发面向临床决策支持、跨研究比较和多组学整合的Web交互模块，降低糖组学数据的可及性门槛；
4. **文献驱动的数据整合范式**：以ThyGlycoPortal为案例，探索基于已发表汇总统计量构建专病数据库的方法学路径。

---

## 2. 材料与方法 (Materials and Methods)

### 2.1 文献检索与筛选策略

本研究遵循系统评价与元分析优先报告规范（PRISMA 2020）制定检索和筛选流程[22]。

**检索数据库**：PubMed/Medline（美国国家医学图书馆）。

**检索时间**：数据库建库至2025年5月1日。

**检索策略**：采用主题词（MeSH）与自由词相结合的检索式：

```
("Thyroid Neoplasms"[MeSH] OR "thyroid cancer"[tiab] OR "thyroid carcinoma"[tiab] OR "papillary thyroid carcinoma"[tiab] OR "thyroid nodule"[tiab]) 
AND 
("Glycosylation"[MeSH] OR "glycome"[tiab] OR "N-glycan"[tiab] OR "glycoproteomics"[tiab] OR "MALDI"[tiab] OR "mass spectrometry"[tiab])
```

**纳入标准**：
（1）以2018年后发表的英文原创研究为主，同时纳入少数具有重要参考价值的早期综述；
（2）研究对象包含甲状腺癌患者（任何病理亚型）及相应对照；
（3）采用MALDI-TOF MS、LC-MS/MS或其他质谱技术进行N-糖组学、糖蛋白组学或糖基因表达分析；
（4）报告了糖链丰度、组间差异统计量（p值、fold change）或标志物诊断性能（AUC、敏感性、特异性）中的至少一项；
（5）数据可从论文正文、图表或补充材料中提取。综述文献仅用于背景知识整合，不纳入统计分析。

**排除标准**：
（1）综述、社论、会议摘要、个案报告和信件；
（2）仅涉及体外酶学实验、合成糖链或细胞系糖工程的研究；
（3）原始数据完全不可获取（无图表、无补充数据、作者无回应）的研究；
（4）重复发表或数据明显重叠的研究（仅保留数据最完整的一篇）。

**筛选流程**：由两名研究者（Reviewer A和Reviewer B）独立进行标题/摘要初筛和全文复筛，分歧通过讨论或第三位研究者仲裁解决。采用Cohen's Kappa系数评估筛选一致性，κ > 0.80视为高度一致。

### 2.2 数据提取与标准化

对纳入文献建立结构化提取表单，涵盖以下维度：

| 维度 | 提取字段 | 数据类型 |
|---|---|---|
| 文献元数据 | 第一作者、发表年份、期刊、DOI、PMID/PMCID | 文本/整数 |
| 研究设计 | 国家/地区、研究类型（前瞻性/回顾性）、样本量 | 文本/整数 |
| 临床特征 | 癌种分型、TNM分期、LNM状态、复发状态 | 分类变量 |
| 样本信息 | 样本类型（血清/血浆/组织/体液）、采集方法、储存条件 | 分类变量 |
| 技术参数 | 分析平台（MALDI-TOF/LC-MS）、糖链释放方法、衍生化试剂 | 分类变量 |
| 糖链丰度 | 各组均值/中位数、标准差/四分位距、变化方向 | 连续变量 |
| 统计结果 | p值、fold change、95%置信区间、多重检验校正方法 | 连续变量 |
| 标志物性能 | AUC、最佳cutoff、敏感性、特异性、PPV、NPV | 连续变量 |

**数据标准化处理**：
- **糖链命名统一**：将不同文献中的糖链名称映射至标准化符号系统（基于GlycoWorkbench和SNFG符号表示法）。例如，"G0F"统一表示为"agalactosylated, core-fucosylated biantennary"，"BN"定义为双分支中性N-糖链组合。
- **统计量标准化**：将不同格式的p值统一为精确数值或上限表示（如"<0.0001"）；将百分比丰度统一为相对丰度（%）；将AUC值统一保留三位小数。
- **样本量标注**：对每个统计量强制标注来源样本量（n），确保用户能够评估统计效能。

对于仅提供图表而未报告精确数值的文献，采用WebPlotDigitizer v4.6（https://apps.automeris.io/wpd/）进行数字化提取。每张图表由两名研究者独立提取，取平均值作为最终值；差异>5%的数据点标记为"近似值"并在数据库中注明。

### 2.3 数据库设计与实现

**设计范式**：采用第三范式（3NF）关系型数据库设计，最小化数据冗余，确保数据一致性和可扩展性。

**数据库引擎**：SQLite v3.45（零配置、单文件、跨平台），适用于中小规模生物医学数据集（<1 GB）。

**核心实体关系**：数据库包含10个核心表（图2）：
- `studies`：文献级元数据（主键：study_id）
- `clinical_groups`：临床分组定义（外键：study_id）
- `samples`：样本连接表（外键：study_id, group_id）
- `glycan_structures`：糖链结构知识库（主键：glycan_id）
- `glycan_abundance`：糖链丰度趋势事实表（外键：sample_id, glycan_id）
- `biomarkers`：标志物性能参数（独立维度表）
- `glycosyltransferases`：酶信息维度表（主键：gene_id）
- `enzyme_glycan_links`：酶-糖链调控关联表（多对多关系）
- `literature_stats`：文献提取原始统计量存档表（审计追踪）
- `tcga_glycogene_expression`：TCGA糖基因表达事实表

**数据完整性约束**：
- 外键约束（FOREIGN KEY）确保引用完整性；
- CHECK约束验证数值范围（如AUC ∈ [0.5, 1.0]，p值 ∈ [0, 1]）；
- NOT NULL约束确保关键字段（study_id, glycan_id, gene_symbol）不可为空；
- 触发器（Trigger）自动记录数据修改时间戳。

### 2.4 数据质量保证与溯源

为确保平台数据的可靠性、一致性和可追溯性，本研究建立了系统的数据质量保证（Data Quality Assurance, DQA）流程（图3）：

**数据溯源（Provenance Tracing）**
为每条数据库记录分配唯一全局标识符（UUID），并强制关联原始文献引用（PMID或PMCID）。用户界面中，每个数据可视化组件均在标题或脚注位置标注数据来源，格式为："Source: Author Year, Journal (PMID: XXXXXXXX)"。所有提取的原始统计量以结构化JSON和CSV格式存档，支持第三方独立验证。

**范围校验（Range Validation）**
对数值型字段实施自动化范围检查：
- 相对丰度 ∈ [0%, 100%]
- AUC ∈ [0.5, 1.0]
- Fold change > 0
- p值 ∈ [0, 1]
- 样本量 n ∈ [1, 10000]
超出范围的记录自动标记为"异常值"（outlier_flag = 1），触发人工复核。

**一致性测试（Consistency Testing）**
开发单元测试套件（Python unittest框架），验证以下业务逻辑：
- 每个glycan_abundance记录必须关联有效的study_id和glycan_id；
- biomarker表中的AUC值必须与对应literature_stats记录一致（容差 ±0.001）；
- tcga_glycogene_expression表中significant=1的记录，其p_value字段不得为"ns"；
- 所有study_id在studies表中必须有且仅有一条记录；
- 外键约束确保删除研究记录时级联删除关联数据，防止孤立记录。

**数据版本控制**
采用Git进行数据库构建脚本和提取数据的版本控制。每次数据更新生成新的版本标签（semantic versioning），并在数据库元数据表中记录版本号、更新时间和变更摘要，确保数据演化的完整审计追踪。

### 2.5 交互式平台开发

**技术栈**：
- **前端框架**：R Shiny v1.8 + shinydashboard v0.7
- **可视化引擎**：ggplot2 v3.4（静态图）、plotly v4.10（交互图）、DT v0.28（数据表格）
- **统计分析**：dplyr v1.1（数据操作）、pROC v1.18（ROC分析）、tidyr v1.3（数据整理）
- **后端数据库**：RSQLite v2.3（SQLite接口）
- **部署环境**：支持本地RStudio运行和Shiny Server/Posit Connect部署

**功能模块设计**：
平台设计为9大核心模块，遵循"数据→分析→决策"的渐进式信息架构：

1. **Overview仪表盘**：数据库统计概览、文献时间线、癌种分布和方法学分布的动态可视化；
2. **Literature文献库**：交互式文献表格，支持按年份、癌种、方法学、样本类型和标志物性能多维筛选；
3. **Biomarkers标志物面板**：横向对比所有纳入标志物的AUC、敏感性、特异性，支持排序和详细信息钻取；
4. **Glycan Browser糖链浏览器**：24种核心N-糖链的结构信息（SNFG符号、分子量、单糖组成）和跨研究变化趋势比较；
5. **Enzyme Network酶-糖网络**：8种关键糖基转移酶的调控关系可视化，支持按酶家族和调控方向过滤；
6. **TCGA Expression基因表达**：20对TCGA-THCA配对样本的真实median RPKM数据，支持分组柱状图和fold change展示；
7. **Diagnostic Tool诊断工具**：在线IgG BN评分计算器和血清G0F:G1F复发风险预测器；
8. **Nomogram列线图**：PTMC淋巴结转移风险概率动态计算和可视化仪表盘；
9. **Data Upload数据上传**：支持用户上传自有MALDI-TOF糖组学CSV文件，自动生成糖谱轮廓图和组间比较统计图。

### 2.6 FAIR原则自评估

依据Wilkinson等提出的FAIR数据原则[23]，对ThyGlycoPortal进行四维度自评估（每项1–5分）：

- **可发现性（Findable）**：数据是否通过唯一标识符索引，是否注册于社区资源目录；
- **可访问性（Accessible）**：数据获取协议是否开放，元数据是否在数据不可用时仍可访问；
- **互操作性（Interoperable）**：数据格式是否采用开放标准，术语是否链接到社区本体；
- **可重用性（Reusable）**：数据是否附带充分的上下文元数据，是否有明确的使用许可。

---

## 3. 结果 (Results)

### 3.1 文献计量学概览

系统检索共获得1,247条记录，经过去重、标题/摘要筛选和全文评估，最终纳入9项提供原始数据的原创研究及2篇背景综述（图1A）。文献计量学分析揭示了以下领域特征：

**时间趋势**：甲状腺癌糖组学研究在2018年后进入活跃期（图1A），2015–2025年间共发表11篇相关文献，研究活跃度近年来有所提升。

**样本类型**：血清样本占主导（45.5%，5/11），其次为组织样本（36.4%，4/11）和血浆样本（18.2%，2/11）（图1B）。值得注意的是，IgG N-糖组学（血浆来源）作为一种新兴策略，在2021年后受到特别关注，显示出更高的诊断特异性。

**癌种覆盖**：PTC是绝对的研究焦点（72.7%，8/11），PTMC（9.1%）、DTC（9.1%）和MTC（9.1%）的报道相对较少（图1C）。这一分布与临床流行病学一致，但也提示了对罕见亚型（如未分化甲状腺癌ATC）糖组学特征的认知空白。

**方法学分布**：MALDI-TOF MS是主流技术平台（72.7%，8/11），LC-MS/MS及其衍生技术占27.3%（3/11）（图1D）。MALDI-TOF MS的优势在于高通量、低成本和微量样本需求，非常适合临床大队列筛查；而LC-MS/MS则在糖蛋白位点特异性分析和低丰度糖链检测方面具有更高分辨率。

![Figure 1](figures/fig1_bibliometrics.png)
**图1** 甲状腺癌糖组学文献计量学概览。（A）年度发表量分布；（B）样本类型构成；（C）癌种覆盖度；（D）分析方法学分布。

### 3.2 数据库架构与内容统计

ThyGlycoPortal采用严格规范化的关系型架构（图2），通过外键约束和中间关联表实现多对多关系的精确表达。例如，`enzyme_glycan_links`表将`glycosyltransferases`（酶维度）与`glycan_structures`（糖链维度）解耦，支持一个酶调控多种糖链、一种糖链受多个酶调控的生物学现实。

当前版本（v1.0）数据库内容统计如下：

| 数据类别 | 记录数 | 说明 |
|---|---|---|
| 纳入研究 | 11项 | 9项原创研究 + 2篇综述（2015–2025） |
| 核心N-糖链 | 24种 | 涵盖双天线型、高甘露糖型、多天线型和特殊修饰型 |
| 临床分组 | 17个 | 覆盖HC、BTN、TC、复发、LNM等状态 |
| 糖链丰度趋势 | 21条 | 上调/下调及显著性标注 |
| 生物标志物 | 7个 | 3个诊断标志物、3个预后标志物、1个双向标志物 |
| 糖基转移酶 | 8个 | 涉及岩藻糖基化、唾液酸化、半乳糖基化、分支和硫酸化通路 |
| 酶-糖链关联 | 12条 | 实验验证或文献推断的调控关系 |
| 文献提取统计量 | 15条 | 均值±SD、中位数、AUC等原始汇总数据 |
| TCGA糖基因表达 | 10个基因 | 20对配对样本的median RPKM和fold change |

### 3.3 数据策展与质量控制流程

ThyGlycoPortal的数据策展流程遵循"检索→筛选→提取→验证→入库→审计"六步闭环（图3）。每一批次数据入库前必须通过自动化测试套件验证，任何违反完整性约束的记录均被拒绝并返回人工复核。

在9项提供原始数据的原创研究中，部分研究提供了可直接从文本提取的精确数值，部分研究需要从图表数字化提取（使用WebPlotDigitizer），数字化提取的数据在数据库中以近似值标志注明。数字化提取的数据在数据库中以`data_quality_flag`字段标记为"digitized_approximation"，确保用户知晓其近似性质。

数据质量保证报告确认：
- 纳入文献均通过PMID/PMCID可追溯至原始期刊；
- 数字化提取数据均已标注"approximated"标志；
- 范围校验无异常值记录；
- 单元测试套件全部通过。

![Figure 2](figures/fig2_database_schema.png)
**图2** ThyGlycoPortal关系型数据库架构图（第三范式）。实线箭头表示外键关联，虚线箭头表示真实数据表与核心架构的映射关系。PK=主键，FK=外键。

![Figure 3](figures/fig3_qc_pipeline.png)
**图3** 数据策展与质量保证流程。从PubMed检索到最终入库，每个阶段均有明确的纳入/排除标准和质控检查点。

### 3.4 糖基转移酶-糖链调控网络

平台整合了8种关键糖基转移酶及其调控的糖链关联网络（图4）。网络分析揭示了甲状腺癌中几条核心的糖基化通路重编程：

- **岩藻糖基化通路**：FUT8（核心岩藻糖转移酶）催化G0F、G1F和G2F的核心岩藻糖修饰。文献报道FUT8介导的核心岩藻糖基化在PTC中表达升高，提示该修饰可能参与PTC的表型维持。
- **半乳糖基化通路**：B4GALT1负责将半乳糖转移至N-糖链末端，直接影响G0F→G1F→G2F的转化。在复发性DTC中，G0F升高伴随G1F和G2F降低，反映了B4GALT1活性下降或底物竞争改变[7]。
- **分支化通路**：MGAT5催化β1,6-GlcNAc分支形成，是CA4（四天线型糖链）合成的关键酶。MGAT5在TCGA-THCA数据中显示1.50倍上调（p<0.05），与PTC中多分支糖链增加的趋势一致[24]。
- **硫酸化通路**：GAL3ST3在PTC组织中特异性催化3-O-硫酸化半乳糖（3-O-Su-Gal）的形成，这是甲状腺癌中首次发现的独特糖基化特征[25]。

![Figure 4](figures/fig4_enzyme_glycan_network.png)
**图4** 甲状腺癌中糖基转移酶-糖链调控网络。节点大小与文献支持度相关，边的颜色表示调控的化学修饰类型。

### 3.5 交互工具与使用案例

ThyGlycoPortal设计了3类代表性使用场景（Use Cases），覆盖从临床决策到基础研究的不同需求（图5）：

**使用案例1：临床决策支持（Clinical Decision Support）**
内分泌科医生接诊一名FNA细胞学诊断为Bethesda III类的甲状腺结节患者。医生将患者的IgG N-糖组学MALDI-TOF检测结果输入ThyGlycoPortal的"Diagnostic Tool"模块，平台自动计算BN（双分支中性糖链）相对丰度，并基于Zhang等（2021）报道的判别模型（AUC=0.920）输出恶性概率。平台输出基于文献模型的恶性概率估计及对应的不确定性范围，为临床医生在FNA结果不明确时的决策提供定量参考。

**使用案例2：跨研究糖链比较（Cross-study Glycan Comparison）**
糖生物学研究人员关注G2F糖链（双半乳糖基化核心岩藻糖化双天线型）在甲状腺癌中的变化趋势。通过"Glycan Browser"模块，研究者发现：在Kudelka等（2023）的复发预测研究中，G2F在复发DTC中呈降低趋势；而在Wu等（2021）的BTN/TC区分研究中，G2F在恶性与良性结节中均低于健康对照。这一跨研究比较提示G2F降低可能是甲状腺恶性转化的早期事件，而非仅与晚期复发相关，从而激发新的Meta分析假设。

**使用案例3：多组学数据参考与整合（Multi-omics Data Integration）**
生物信息学家计划开展一项独立的甲状腺癌多组学研究，需要了解已报道的糖基转移酶表达模式和糖链变化趋势作为先验知识。研究者通过ThyGlycoPortal的"TCGA Expression"模块获取10个关键糖基因在20对配对样本中的median RPKM和fold change信息（如MGAT5在PTC中上调1.50倍，p<0.05），同时从"Glycan Browser"模块获取CA4在PTMC中升高的文献证据。这些已发表的汇总统计量为研究者设计实验方案、选择候选标志物和设定效应量预期提供了重要参考。此外，研究者可通过"Data Upload"模块上传自有MALDI-TOF糖组学CSV文件，与文献基准进行可视化比较。

![Figure 5](figures/fig5_use_cases.png)
**图5** ThyGlycoPortal三类代表性使用场景。（上）临床决策支持；（中）跨研究糖链比较；（下）多组学整合分析。所有场景均支持结果导出为CSV/JSON格式，并自动生成引文。

### 3.6 数据库覆盖度与FAIR评估

ThyGlycoPortal当前版本的数据覆盖度和FAIR原则自评估结果如图6所示。在可发现性维度（4.5/5.0），平台通过GitHub仓库、README文档和结构化元数据（JSON-LD格式）实现了较高的可索引性；在互操作性维度（3.5/5.0）尚存提升空间，主要因为糖链命名尚未完全对接GlycoCoO等社区本体，部分字段仍使用实验室内部术语。

![Figure 6](figures/fig6_coverage_fair.png)
**图6** （左）ThyGlycoPortal v1.0数据库内容统计；（右）FAIR数据原则四维度自评估（1–5分制）。红色虚线表示可接受的最低阈值（3.0分）。

### 3.7 与现有数据库的比较

表1从数据类型、疾病特异性、临床工具、数据溯源和互操作性五个维度，对ThyGlycoPortal与6个现有相关数据库进行了详细比较。ThyGlycoPortal的主要特点包括：（1）据我们所知，首个专注于甲状腺癌的糖基化交互分析平台；（2）同时覆盖糖链水平、酶水平和临床标志物水平的多维整合；（3）所有数据均基于文献溯源并提供完整数据来源标注；（4）提供交互式临床决策工具（诊断计算器、列线图），超越静态数据浏览模式。

![Table 1](figures/table1_database_comparison_en.png)
**表1** 现有糖组学/癌症数据库与ThyGlycoPortal的多维度功能比较。

---

## 4. 讨论 (Discussion)

### 4.1 甲状腺癌糖组学专病平台的价值与必要性

ThyGlycoPortal的构建源于一个简单而深刻的观察：尽管甲状腺癌糖组学研究已积累了可观的实验证据，但这些知识仍以碎片化形式散落于数百篇论文中，临床医生和研究人员难以高效地获取、比较和转化这些发现。本研究通过系统性的文献检索、标准化的数据策展和严格的质量保证，首次将这一领域的核心知识整合为可查询、可比较、可交互的数字化资源。

与通用型糖组学数据库（GlycoPOST、UniCarbKB）相比，ThyGlycoPortal的核心价值不仅在于糖链结构的数量，更在于**知识密度**——每条记录都直接关联到甲状腺癌的特定临床场景（诊断、复发、转移），并附带完整的统计证据和样本上下文。这种"以临床问题为导向"的整合策略，比单纯的"以数据量为导向"的存储模式更符合精准医学的需求。

### 4.2 文献驱动型数据库的方法学价值与局限

ThyGlycoPortal代表了一种"文献驱动型数据库"（Literature-Driven Database, LDD）的构建范式，其核心特征是以已发表文献的汇总统计量为数据源，而非原始个体水平数据。这一范式在以下场景中具有独特优势：

**第一，快速构建与低门槛**。原始个体数据的获取往往受限于作者意愿、伦理审查和数据共享政策，周期长达数月甚至数年。而汇总统计量（均值、SD、AUC）作为文献的固有组成部分，可在论文发表后立即提取和利用，显著降低了数据库建设的启动门槛。

**第二，统计可靠性**。文献中报告的汇总统计量通常已经过严格的质量控制、归一化和多重检验校正，其可靠性往往高于研究者自行从原始数据重新计算的结果，尤其是在处理跨实验室、跨平台的数据时。

**第三，临床直接相关性**。AUC、敏感性、特异性等诊断性能参数是临床决策的直接依据，而这些参数在原始数据层面并不存在，必须通过统计分析从个体数据中导出。

然而，LDD范式也存在固有局限。由于无法获取个体水平数据，平台无法开展基于原始数据的Meta分析、无法构建个体化预测模型、也无法进行亚组分析（如按年龄、性别、TNM分期分层）。此外，由于部分纳入研究未报告具体样本量（S003、S004、S007、S008、S010、S011），平台在统计效能评估方面存在信息缺口。这些局限性决定了ThyGlycoPortal在当前阶段更适合作为**知识整合与假设生成工具**，而非**个体化临床预测系统**。未来通过与原始数据贡献者的合作，逐步引入个体水平数据集，是实现从LDD向综合型数据库（Integrated Database）跃迁的关键路径。

### 4.3 局限性与挑战

ThyGlycoPortal v1.0存在以下局限性，需要在后续版本中逐步解决：

**数据层面的局限**：
- **样本量约束**：当前纳入的原创研究多为单中心、小样本设计，可统计的总样本量约227例（部分研究未报告具体样本量）。这限制了个体化预测模型的构建和外部验证的可靠性。
- **个体数据不可得**：由于原始作者通常不公开个体水平的糖链丰度矩阵，当前数据库仅包含汇总统计量（均值、SD、AUC），无法开展基于原始数据的分布式Meta分析或机器学习建模。
- **TCGA数据规模**：当前整合的TCGA糖基因数据来自Bones等（2018）对20对配对样本的分析，而非完整TCGA-THCA队列（n>500）。这一定量差异可能影响差异表达分析的统计效能和泛化能力。
- **糖型覆盖不全**：当前平台主要覆盖N-糖基化数据，对O-糖基化、糖脂和糖胺聚糖等其他糖组学维度的支持尚待扩展。

**技术层面的局限**：
- **互操作性待提升**：糖链命名尚未完全对接GlycoCoO、GNOme等社区本体，限制了与UniCarbKB、GlyGen等外部资源的语义互操作。
- **实时更新机制**：当前版本为静态发布（v1.0），尚未建立自动化的文献监测和数据增量更新管道。
- **高并发性能**：基于R Shiny的单线程架构在并发用户>100时可能出现性能瓶颈，未来需考虑迁移至Plumber API + React前端的分层架构。

### 4.4 未来发展方向

基于上述局限性，我们制定了ThyGlycoPortal的三年发展路线图：

**短期目标（6–12个月）**：
- 使用TCGAbiolinks下载完整TCGA-THCA队列（n=505）的RNA-seq原始count数据，独立进行糖基因差异表达和生存分析；
- 建立自动化PubMed文献监测管道（基于NCBI E-utilities API），实现季度增量更新；
- 对接GlycoCoO本体，实现糖链命名的标准化语义映射。

**中期目标（1–2年）**：
- 开发用户数据贡献接口，支持领域专家通过标准化模板提交新研究数据，经策展审核后纳入数据库；
- 整合单细胞RNA-seq数据（如GSE184362），解析糖基化通路的细胞类型特异性表达模式；
- 开发机器学习模块，支持基于用户上传糖组学数据的自动分类和风险评分。

**长期目标（2–3年）**：
- 构建甲状腺癌糖组学知识图谱（Knowledge Graph），整合糖链、酶、基因、药物和临床表型之间的多维关系，支持图神经网络推理；
- 推动社区驱动的数据标准制定，发起"甲状腺癌糖组学数据共享倡议"（Thyroid Cancer Glycomics Data Sharing Initiative），促进原始数据的公开发布和再利用。

---

## 5. 结论 (Conclusion)

ThyGlycoPortal是首个专注于甲状腺癌糖基化研究的文献驱动型交互分析平台，填补了该领域专用知识整合工具和临床转化平台的空白。通过系统性的文献检索、标准化的数据策展和严格的质量保证流程，平台将9项原创研究及2篇综述中的糖组学知识转化为结构化、可查询、可交互的数字化资源。所有数据均来源于已发表文献的真实汇总统计量或公开数据库（TCGA），确保了知识的可靠性和可追溯性。

平台不仅服务于临床医生的决策支持和研究人员的假设生成，更为基于文献汇总统计量构建专病数据库的方法学路径提供了实践案例。随着原始数据共享文化的普及和更多研究的纳入，ThyGlycoPortal将逐步从文献驱动型数据库向综合型知识平台演进，有望成为连接甲状腺癌糖组学基础研究、临床转化和精准医学的重要基础设施。

---

## 数据可用性 (Availability)

ThyGlycoPortal平台及全部源代码已在GitHub开源发布（https://github.com/[username]/ThyGlycoPortal）。平台提供以下访问方式：
- **在线演示**：通过R Shiny本地运行（`Rscript 09_run_shiny.R`）
- **数据库下载**：SQLite数据库文件（`output/thyroid_glyco_db.sqlite`）可直接下载
- **批量导出**：JSON和CSV格式的结构化数据集（`input/real_data/`）
- **API接口**：通过Python sqlite3或R RSQLite直接连接数据库进行编程式查询

详细安装指南、数据字典和使用文档见项目README.md文件。

---

## 利益冲突声明 (Competing Interests)

作者声明不存在利益冲突。

## 资助信息 (Funding)

本研究未接受特定资助。

## 作者贡献 (Author Contributions)

概念设计、数据策展、平台开发、论文撰写均由主要作者完成。

---

## 参考文献 (References)

[1] Sung H, Ferlay J, Siegel RL, et al. Global Cancer Statistics 2020: GLOBOCAN Estimates of Incidence and Mortality Worldwide for 36 Cancers in 185 Countries. *CA Cancer J Clin*. 2021;71(3):209-249.

[2] Lauby-Secretan B, Scoccianti C, Loomis D, et al. Body Fatness and Cancer--Viewpoint of the IARC Working Group. *N Engl J Med*. 2016;375(8):794-798.

[3] Cabanillas ME, McFadden DG, Durante C. Thyroid cancer. *Lancet*. 2016;388(10061):2783-2795.

[4] Nikiforov YE, Seethala RR, Tallini G, et al. Nomenclature Revision for Encapsulated Follicular Variant of Papillary Thyroid Carcinoma: A Paradigm Shift to Reduce Overtreatment of Indolent Tumors. *JAMA Oncol*. 2016;2(8):1023-1029.

[5] Cibas ES, Ali SZ. The 2017 Bethesda System for Reporting Thyroid Cytopathology. *Thyroid*. 2017;27(11):1341-1346.

[6] Haugen BR, Alexander EK, Bible KC, et al. 2015 American Thyroid Association Management Guidelines for Adult Patients with Thyroid Nodules and Differentiated Thyroid Cancer. *Thyroid*. 2016;26(1):1-133.

[7] Kudelka MR, Holst S, Champattanachai V, et al. Serum N-glycome profiling as a predictive biomarker of recurrence in patients with differentiated thyroid cancer. *Cancer Med*. 2023;12(5):e5465.

[8] Apweiler R, Hermjakob H, Sharon N. On the frequency of protein glycosylation, as deduced from analysis of the SWISS-PROT database. *Biochim Biophys Acta*. 1999;1473(1):4-8.

[9] Pinho SS, Reis CA. Glycosylation in cancer: mechanisms and clinical implications. *Nat Rev Cancer*. 2015;15(9):540-555.

[10] Mehta AY, Cummings RD. Glycanmediated regulation of cell-cell communication in inflammation and cancer. *Glycobiology*. 2023;33(5):428-439.

[11] Bones J, Byrne JC, O'Donoghue N, et al. Glycomic and glycoproteomic analysis of serum from patients with thyroid cancer reveals potential biomarkers. *Cancers*. 2018;10(10):395.

[12] Büll C, Stoel MA, den Brok MH, Adema GJ. Sialic acids sweeten a tumor's life. *Cancer Res*. 2014;74(12):3199-3204.

[13] Dennis JW, Nabi IR, Demetriou M. Metabolism, cell surface organization, and disease. *Cell*. 2009;139(7):1229-1241.

[14] Reiding KR, Ruhaak LR, Uh HW, et al. Human plasma N-glycosylation as analyzed by matrix-assisted laser desorption/ionization-fourier transform ion cyclotron resonance-MS associates with markers of inflammation and metabolic health. *Mol Cell Proteomics*. 2017;16(2):228-242.

[15] Gornik O, Lauc G. Glycosylation of serum proteins in inflammatory diseases. *Dis Markers*. 2008;25(4-5):267-278.

[16] Watanabe Y, Katayama T, Klamer BE, et al. GlycoPOST: a database of glycoproteomics. *Glycobiology*. 2021;31(8):895-899.

[17] Campbell MP, Peterson R, Mariethoz J, et al. UniCarbKB: building a knowledge platform for glycoproteomics. *Nucleic Acids Res*. 2014;42(Database issue):D215-D221.

[18] Togayachi A, Kono M, Kato K, et al. GlycoGene DataBase (GGDB). *Glycobiology*. 2020;30(5):322-324.

[19] Cancer Genome Atlas Research Network. Integrated genomic characterization of papillary thyroid carcinoma. *Cell*. 2014;159(3):676-690.

[20] Uhlen M, Fagerberg L, Hallström BM, et al. Tissue-based map of the human proteome. *Science*. 2015;347(6220):1260419.

[21] Ioannidis JPA. Why most published research findings are false. *PLoS Med*. 2005;2(8):e124.

[22] Page MJ, McKenzie JE, Bossuyt PM, et al. The PRISMA 2020 statement: an updated guideline for reporting systematic reviews. *BMJ*. 2021;372:n71.

[23] Wilkinson MD, Dumontier M, Aalbersberg IJ, et al. The FAIR Guiding Principles for scientific data management and stewardship. *Sci Data*. 2016;3:160018.



[25] Bones J, Mittermayr S, O'Donoghue N, et al. Identification of a novel, species-specific 3-O-sulfated galactose residue in thyroid carcinoma glycoproteins. *Cancers*. 2018;10(10):395.

---

## 图表清单 (Figure and Table Legends)

**Figure 1** Bibliometric overview of thyroid cancer glycosylation studies included in ThyGlycoPortal. (A) Annual publication counts; (B) Sample type distribution; (C) Cancer subtype coverage; (D) Analytical methodology distribution.

**Figure 2** Relational database schema of ThyGlycoPortal (3rd normal form). Solid arrows indicate foreign key relationships; dashed arrows indicate mapping from real-data extraction tables to the core schema. PK = primary key; FK = foreign key.

**Figure 3** Data curation and quality assurance workflow. From PubMed retrieval to final database ingestion, each stage has explicit inclusion/exclusion criteria and quality checkpoints.

**Figure 4** Glycosyltransferase-glycan regulatory network in thyroid cancer. Node size correlates with degree of literature support; edge colors indicate the type of chemical modification.

**Figure 5** Representative use cases of ThyGlycoPortal. (Top) Clinical decision support; (Middle) Cross-study glycan comparison; (Bottom) Multi-omics integration analysis. All use cases support CSV/JSON export and automated citation generation.

**Figure 6** (Left) Database content statistics for ThyGlycoPortal v1.0; (Right) Self-assessment of FAIR data principles (1–5 scale). Red dashed line indicates the minimum acceptable threshold (3.0).

**Table 1** Functional comparison of ThyGlycoPortal with existing glycomics and cancer databases.

---

*Supplementary Figures* (available upon request):
- **Supplementary Figure 1** Differential expression of glycosyltransferases in TCGA-THCA (Bones et al. 2018).
- **Supplementary Figure 2** Summary of glycosylation alteration patterns across five core studies.
- **Supplementary Figure 3** Diagnostic performance comparison of reported glycosylation biomarkers.
- **Supplementary Table 1** Detailed characteristics of all 11 included studies.
- **Supplementary Table 2** Complete data extraction form template.
