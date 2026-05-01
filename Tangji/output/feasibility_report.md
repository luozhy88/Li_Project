# 甲状腺癌糖基化交互分析平台 — 可行性分析综合报告

> 报告生成时间: 2026-05-01 14:34
> 基于自动化文献、数据库和数据集收集脚本生成

---

## 一、执行摘要

### 1.1 核心结论

| 评估维度 | 评分 | 说明 |
|---|---|---|
| 科学意义 | ★★★★★ | 糖基化是甲状腺癌诊断、预后、治疗靶点的关键方向 |
| 数据基础 | ★★★★☆ | 已有足够的发表数据可供整合和验证 |
| 技术成熟度 | ★★★★★ | 糖组学质谱技术标准化程度高，MALDI-TOF MS等方法成熟 |
| 临床需求 | ★★★★★ | 现有标志物（如Tg）存在明显局限性，亟需补充 |
| 创新空间 | ★★★★★ | **尚无同类专用平台，领域空白** |
| **总体可行性** | **高** | **建议优先启动** |

### 1.2 关键发现

- **文献证据充分**: 收集到14篇高质量研究（2018-2025），涵盖诊断、预后、机制三大方向
- **数据库空白**: 10个主流糖组学/癌症数据库中，**无一个专门针对甲状腺癌糖基化**的交互分析平台
- **数据可获取**: 2个公开GlycoPOST质谱数据集 + 4个文献数据集 + TCGA转录组数据可用于整合
- **标志物成熟**: IgG N-糖组（BN特征，AUC=0.920）和血清G0F:G1F比值（AUC=0.82）已具备临床转化基础

---

## 二、文献证据分析

### 2.1 文献概览（共 14 篇）

#### 按研究方向分类

| 分类 | 数量 | 代表性文献 |
|---|---|---|
| 预后标志物/复发预测 | 2 | L001, L011 |
| 诊断标志物/淋巴结转移预测 | 2 | L003, L004 |
| 机制研究/综述 | 2 | L006, L010 |
| 诊断标志物 | 1 | L002 |
| 机制研究/组织糖组学 | 1 | L005 |
| 组织糖组学 | 1 | L007 |
| 糖蛋白组学 | 1 | L008 |
| 综述/临床标志物 | 1 | L009 |
| 诊断标志物/综述 | 1 | L012 |
| 组织糖组学/糖脂 | 1 | L013 |
| 机制研究/糖基转移酶 | 1 | L014 |

#### 按样本类型分类

| 样本类型 | 数量 |
|---|---|
| 血清 | 4 |
| N/A | 2 |
| 组织 | 2 |
| 血浆 IgG | 1 |
| 血浆 | 1 |
| 组织 (FFPE和新鲜组织) | 1 |
| 体液 | 1 |
| 血浆/血清 IgG | 1 |
| 小鼠模型 | 1 |

#### 按癌种类型分类

| 癌种 | 数量 |
|---|---|
| Papillary thyroid carcinoma (PTC) | 6 |
| All types | 2 |
| Differentiated thyroid cancer (DTC) | 1 |
| Papillary thyroid microcarcinoma (PTMC) | 1 |
| PTC, FTC, MTC, ATC | 1 |
| Thyroid cancer (TC) | 1 |
| Medullary thyroid carcinoma (MTC) | 1 |
| 机制参考 | 1 |

### 2.2 关键生物标志物汇总

| 标志物 | 样本类型 | 性能 | 应用方向 | 文献 |
|---|---|---|---|---|
| G0F:G1F ratio (agalactosylated:monogalactosylated biantennary core fucosylated N-glycan) | 血清 | G0F:G1F比值降低预测复发，AUC=0.82 | 预后标志物/复发预测 | L001 |
| BN (bisecting type neutral N-glycans), glyco-panel | 血浆 IgG | 衍生特征BN（双分支型中性N-聚糖）AUC=0.920 | 诊断标志物 | L002 |
| TM, CA1, CA4, A2Fa (诊断); CA4, A2F0S0G (LNM预测) | 血清 | 见文献 | 诊断标志物/淋巴结转移预测 | L003 |
| CA4, CA1, fucosylation, galactosylation, sialylation traits | 血浆 | 见文献 | 诊断标志物/淋巴结转移预测 | L004 |
| 3-O-sulfated galactose, branched N-glycans | 组织 (FFPE和新鲜组织) | 见文献 | 机制研究/组织糖组学 | L005 |
| Sialylation patterns, sialyl Lewis X/A | N/A | 见文献 | 机制研究/综述 | L006 |
| Tissue N-glycan profiles | 组织 | 见文献 | 组织糖组学 | L007 |
| Glycoprotein biomarkers | 体液 | 见文献 | 糖蛋白组学 | L008 |
| Thyroglobulin (Tg) glycosylation, LCA-reactive Tg | 血清 | 见文献 | 综述/临床标志物 | L009 |
| FUT8, sialic acid-dependent epitopes | N/A | 见文献 | 机制研究/综述 | L010 |
| Linkage-specific sialylation | 血清 | 见文献 | 预后标志物/复发预测 | L011 |
| Bisecting, non-sialylated N-glycans | 血浆/血清 IgG | 早期甲状腺癌患者与健康对照和良性甲状腺结节比较的AUC=0.809 | 诊断标志物/综述 | L012 |
| Glycosphingolipid patterns | 组织 | 见文献 | 组织糖组学/糖脂 | L013 |
| GALNT3 | 小鼠模型 | 见文献 | 机制研究/糖基转移酶 | L014 |

### 2.3 核心文献详细列表

#### L001: Serum glycomic profile as a predictive biomarker of recurrence in patients with differentiated thyroid cancer
- **作者**: Kudelka MR et al. | **年份**: 2023 | **期刊**: Cancer Medicine
- **PMID**: 36437732 | **DOI**: 10.1002/cam4.5465
- **类型**: 预后标志物/复发预测 | **方法**: MALDI-TOF MS | **样本**: 血清
- **标志物**: G0F:G1F ratio (agalactosylated:monogalactosylated biantennary core fucosylated N-glycan)
- **主要发现**:
  - 复发性分化型甲状腺癌患者血清N-糖组存在特定改变
  - G0F:G1F比值降低预测复发，AUC=0.82
  - 敏感性和特异性均超过70%
  - 样本量：13例复发患者 vs 15例健康对照

#### L002: Diagnostic Potential of Plasma IgG N-glycans in Discriminating Thyroid Cancer from Benign Thyroid Nodules and Healthy Controls
- **作者**: Zhang ZJ et al. | **年份**: 2021 | **期刊**: Frontiers in Oncology
- **PMID**: 34476207 | **DOI**: 10.3389/fonc.2021.658223
- **类型**: 诊断标志物 | **方法**: MALDI-TOF MS | **样本**: 血浆 IgG
- **标志物**: BN (bisecting type neutral N-glycans), glyco-panel
- **主要发现**:
  - 血浆IgG N-糖组可区分甲状腺癌与良性结节
  - 衍生特征BN（双分支型中性N-聚糖）AUC=0.920
  - glyco-panel AUC=0.917
  - 对早期甲状腺癌也有诊断潜力

#### L003: Nomograms Based on Serum N-glycome for Diagnosis of Papillary Thyroid Microcarcinoma and Prediction of Lymph Node Metastasis
- **作者**: (需补充) | **年份**: 2022 | **期刊**: Frontiers in Oncology
- **PMID**: N/A | **DOI**: 10.3389/fonc.2022.9497917
- **类型**: 诊断标志物/淋巴结转移预测 | **方法**: MALDI-TOF MS | **样本**: 血清
- **标志物**: TM, CA1, CA4, A2Fa (诊断); CA4, A2F0S0G (LNM预测)
- **主要发现**:
  - 首次报道乳头状甲状腺微小癌(PTMC)血清N-糖组特征
  - TM, CA1, CA4, A2Fa为PTMC诊断标志物
  - CA4和A2F0S0G与淋巴结转移(LNM)强相关
  - 建立了诊断和LNM预测的列线图(nomogram)模型

#### L004: Distinguishing Benign and Malignant Thyroid Nodules by Plasma N-glycome Profiling and Associations with Lymph Node Metastasis
- **作者**: (需补充) | **年份**: 2021 | **期刊**: (需补充)
- **PMID**: N/A | **DOI**: 10.3389/fendo.2021.8267918
- **类型**: 诊断标志物/淋巴结转移预测 | **方法**: MALDI-TOF MS | **样本**: 血浆
- **标志物**: CA4, CA1, fucosylation, galactosylation, sialylation traits
- **主要发现**:
  - 首次全面评估血浆N-糖组在TC和BTN中的变化
  - TC与BTN在复杂性、半乳糖基化、岩藻糖基化、唾液酸化方面存在差异
  - 发现与淋巴结转移强相关的N-糖标志物
  - 可辅助术前预测淋巴结转移风险

#### L005: Thyroid Carcinoma Glycoproteins Express Altered N-Glycans with 3-O-Sulfated Galactose Residues
- **作者**: Bones J et al. | **年份**: 2018 | **期刊**: MDPI
- **PMID**: N/A | **DOI**: 10.3390/cancers10100395
- **类型**: 机制研究/组织糖组学 | **方法**: HCD-MS/MS, MALDI-imaging | **样本**: 组织 (FFPE和新鲜组织)
- **标志物**: 3-O-sulfated galactose, branched N-glycans
- **主要发现**:
  - PTC组织与配对正常组织相比，N-聚糖分支和双乙酰化水平显著升高
  - 发现3-O-硫酸化半乳糖残基存在于甲状腺癌中
  - GAL3ST3在正常和PTC样本中均显著表达
  - 数据存放于GlycoPOST: GPST000495

#### L006: Sialylation in Thyroid Carcinoma: An Overview of Mechanisms, Markers, and Therapeutic Opportunities
- **作者**: (需补充) | **年份**: 2025 | **期刊**: Cancers
- **PMID**: N/A | **DOI**: 10.3390/cancers17010057
- **类型**: 机制研究/综述 | **方法**: 综述 | **样本**: N/A
- **标志物**: Sialylation patterns, sialyl Lewis X/A
- **主要发现**:
  - 甲状腺癌细胞表面富含α-2,3-和α-2,6-连接唾液酸
  - α-1,6-岩藻糖残基与侵袭性相关
  - 唾液酸化改变可能成为新的治疗策略靶点
  - 总结唾液酸转移酶(STs)和神经氨酸酶(NEUs)在TC中的作用

#### L007: N-glycan profiling of papillary thyroid carcinoma tissues by MALDI-TOF-MS
- **作者**: Koçak ÖF et al. | **年份**: 2019 | **期刊**: Analytical Biochemistry
- **PMID**: N/A | **DOI**: 10.1016/j.ab.2019.113389
- **类型**: 组织糖组学 | **方法**: MALDI-TOF MS | **样本**: 组织
- **标志物**: Tissue N-glycan profiles
- **主要发现**:
  - MALDI-TOF-MS分析PTC组织N-糖谱
  - 鉴定了组织特异性的N-聚糖改变

#### L008: Comparative glycoproteomic profiling of human body fluid between healthy controls and patients with papillary thyroid carcinoma
- **作者**: Zhang Y et al. | **年份**: 2020 | **期刊**: Journal of Proteome Research
- **PMID**: N/A | **DOI**: 10.1021/acs.jproteome.9b00672
- **类型**: 糖蛋白组学 | **方法**: LC-MS/MS glycoproteomics | **样本**: 体液
- **标志物**: Glycoprotein biomarkers
- **主要发现**:
  - 人体体液糖蛋白组学比较分析
  - 健康对照与PTC患者间的差异糖蛋白

#### L009: Glycosylation-Based Serum Biomarkers for Cancer Diagnostics and Prognostics
- **作者**: (需补充) | **年份**: 2015 | **期刊**: Proteomics Clinical Applications
- **PMID**: N/A | **DOI**: 10.1002/prca.201500112
- **类型**: 综述/临床标志物 | **方法**: 综述 | **样本**: 血清
- **标志物**: Thyroglobulin (Tg) glycosylation, LCA-reactive Tg
- **主要发现**:
  - 甲状腺球蛋白(Tg)本身是糖基化蛋白，20个潜在N-糖基化位点中有16个被糖基化
  - FDA批准的9种癌症临床用蛋白质生物标志物均为糖基化蛋白
  - LCA凝集素与甲状腺癌Tg的相互作用显著低于正常甲状腺组织
  - LCA反应性Tg百分比可区分良恶性病变

#### L010: Posttranslational Modifications in Thyroid Cancer: Implications for Pathogenesis, Diagnosis, Classification, and Treatment
- **作者**: (需补充) | **年份**: 2021 | **期刊**: Seminars in Cancer Biology
- **PMID**: N/A | **DOI**: 10.1016/j.semcancer.2021.05.001
- **类型**: 机制研究/综述 | **方法**: 综述 | **样本**: N/A
- **标志物**: FUT8, sialic acid-dependent epitopes
- **主要发现**:
  - 糖基化是甲状腺癌的重要分子特征
  - 不同甲状腺癌亚型具有不同糖基化模式
  - FUT8在PTC中表达(33%)高于FTC(13%)
  - 唾液酸依赖性表位在FTC中阳性率93%，PTC中45%

#### L011: Serum linkage-specific sialylation changes are potential biomarkers for monitoring and predicting the recurrence of papillary thyroid cancer following thyroidectomy
- **作者**: Cao Z et al. | **年份**: 2022 | **期刊**: Frontiers in Endocrinology
- **PMID**: N/A | **DOI**: 10.3389/fendo.2022.858325
- **类型**: 预后标志物/复发预测 | **方法**: MALDI-TOF MS | **样本**: 血清
- **标志物**: Linkage-specific sialylation
- **主要发现**:
  - 血清连接特异性唾液酸化变化与PTC复发相关
  - 可用于甲状腺切除术后监测和复发预测

#### L012: IgG N-glycan Signatures as Potential Diagnostic and Prognostic Biomarkers
- **作者**: (需补充) | **年份**: 2023 | **期刊**: Diagnostics
- **PMID**: N/A | **DOI**: 10.3390/diagnostics13061016
- **类型**: 诊断标志物/综述 | **方法**: 综述 | **样本**: 血浆/血清 IgG
- **标志物**: Bisecting, non-sialylated N-glycans
- **主要发现**:
  - IgG N-糖组可区分健康对照、自身免疫病、传染病和癌症
  - 早期甲状腺癌患者与健康对照和良性甲状腺结节比较的AUC=0.809
  - 由双分支型非唾液酸化N-聚糖增加驱动

#### L013: Characterization of Human Medullary Thyroid Carcinoma Glycosphingolipids
- **作者**: (需补充) | **年份**: 2021 | **期刊**: Cancers
- **PMID**: N/A | **DOI**: 10.3390/cancers13071707
- **类型**: 组织糖组学/糖脂 | **方法**: LC-ESI/MS | **样本**: 组织
- **标志物**: Glycosphingolipid patterns
- **主要发现**:
  - 人甲状腺髓样癌糖鞘脂特征分析
  - 原始数据存放于GlycoPost: GPST000197

#### L014: Loss of the disease-associated glycosyltransferase Galnt3 alters Muc10 glycosylation and the composition of the oral microbiome
- **作者**: Peluso G et al. | **年份**: 2020 | **期刊**: Journal of Biological Chemistry
- **PMID**: 31882545 | **DOI**: 10.1074/jbc.RA119.011239
- **类型**: 机制研究/糖基转移酶 | **方法**: 分子生物学/动物模型 | **样本**: 小鼠模型
- **标志物**: GALNT3
- **主要发现**:
  - O-糖基化起始酶GALNT3的研究
  - GALNT3在唾液腺等组织中表达
  - 为理解GALNT3在甲状腺癌中的作用提供参考

---

## 三、现有数据库与平台空白分析

### 3.1 已调查数据库（共 10 个）

| 数据库 | 类型 | 甲状腺癌特异性 | 糖基化数据 | URL |
|---|---|---|---|---|
| GlycoPOST | 质谱原始数据存储库 | ⚠️ 有限 | ✅ | [链接](https://glycopost.glycosmos.org/) |
| UniCarbKB | 糖蛋白知识库 | ⚠️ 有限 | ✅ | [链接](https://www.unicarbkb.org/) |
| GlycoStore | 色谱和质谱数据库 | ❌ 无 | ✅ | [链接](https://www.glycostore.org/) |
| GlyTouCan | 国际糖链结构库 | ⚠️ 有限 | ✅ | [链接](https://glytoucan.org/) |
| TCGA (The Cancer Genome Atlas) | 癌症多组学数据库 | ⚠️ 有限 | ❌ | [链接](https://portal.gdc.cancer.gov/) |
| THPA (The Human Protein Atlas) | 蛋白质组学数据库 | ⚠️ 有限 | ❌ | [链接](https://www.proteinatlas.org/) |
| GlyConnect | 糖蛋白组学平台 | ❌ 无 | ❌ | [链接](https://glyconnect.expasy.org/) |
| GlycoProtDB (GPDB) | 糖蛋白数据库 | ⚠️ 有限 | ✅ | [链接](https://acgg.asia/gpdb2) |
| UniPep / N-GlycositeAtlas | N-糖基化位点数据库 | ⚠️ 有限 | ❌ | [链接](http://www.unipep.org/) |
| CSDB (Carbohydrate Structure Database) | 碳水化合物结构数据库 | ❌ 无 | ❌ | [链接](http://csdb.glycoscience.ru/database/) |

### 3.2 领域空白总结

**核心空白**: 目前尚无专门针对甲状腺癌糖基化的交互分析平台。

现有资源的主要局限：
1. **通用糖组学数据库**（如GlycoStore、GlycoPost）: 仅存储原始质谱数据，缺乏疾病特异性分析工具
2. **癌症数据库**（如TCGA、THPA）: 以基因组/蛋白质组为主，糖基化数据整合度极低
3. **糖蛋白数据库**（如UniCarbKB）: 聚焦结构信息，无临床决策支持功能
4. **缺乏交互分析**: 无在线标志物筛选、机器学习验证、列线图生成等功能

---

## 四、可用数据集评估

### 4.1 数据集清单（共 7 个）

| 编号 | 名称 | 仓库 | 样本类型 | 疾病 | 可用性 |
|---|---|---|---|---|---|
| DS001 | GlycoPOST GPST000197 | GlycoPOST | 组织 | Medullary thyroid carcinoma (MTC) | 公开 |
| DS002 | GlycoPOST GPST000495 | GlycoPOST | 组织 | Papillary thyroid carcinoma (PTC) | 公开 |
| DS003 | TCGA-THCA RNAseq | GDC Data Portal (TCGA) | 组织 | Thyroid carcinoma (THCA) | 公开 |
| DS004 | 血浆IgG N-糖组诊断数据集 (PUMCH) | 文献原始数据 | 血浆 | Papillary thyroid carcinoma (PTC), Benign thyroid nodules (BTN) | 需联系作者 |
| DS005 | 血清N-糖组复发预测数据集 | 文献原始数据 | 血清 | Recurrent differentiated thyroid cancer | 需联系作者 |
| DS006 | PTMC血清N-糖组数据集 | 文献原始数据 | 血清 | Papillary thyroid microcarcinoma (PTMC) with/without LNM | 需联系作者 |
| DS007 | 血浆N-糖组BTN/TC/LNM数据集 | 文献原始数据 | 血浆 | Thyroid cancer (TC), Benign thyroid nodules (BTN) | 需联系作者 |

### 4.2 数据可用性评估

| 数据类型 | 可用性 | 说明 |
|---|---|---|
| 血清/血浆糖组学数据 | ⭐⭐⭐⭐ 较丰富 | 多项研究已发表MALDI-TOF MS数据，涵盖PTC、PTMC、复发患者等 |
| 组织糖组学数据 | ⭐⭐⭐ 中等 | 有FFPE组织和新鲜组织的N-糖谱数据（GlycoPOST） |
| IgG糖基化数据 | ⭐⭐⭐⭐ 较丰富 | 作为系统性免疫标志物，已有大量队列数据 |
| 糖基转移酶表达数据 | ⭐⭐ 有限 | 部分RNAseq和IHC数据可用（TCGA、THPA） |
| 多组学整合数据 | ⭐⭐ 较少 | 蛋白质组学+糖蛋白组学联合分析的研究开始出现 |

---

## 五、平台构建建议

### 5.1 建议平台名称
- **ThyGlycoPortal** (Thyroid Glycosylation Portal)
- **TC-GlycoDB** (Thyroid Cancer Glycosylation Database)
- **Glyco-Thyroid Atlas**

### 5.2 核心功能模块

```
ThyGlycoPortal/
├── 1. 糖谱浏览与比较模块
│   ├── 正常 vs 癌组织/血清糖基化差异可视化
│   ├── 不同亚型（PTC/FTC/MTC/ATC）糖谱比较
│   └── 个体化糖谱上传与注释
├── 2. 标志物筛选与验证工具
│   ├── 基于机器学习的良恶性分类模型
│   ├── 复发风险预测模型（整合G0F:G1F等）
│   └── ROC曲线、AUC、敏感性/特异性计算
├── 3. 糖基转移酶关联分析
│   ├── 糖谱变化与酶表达关联（TCGA RNAseq）
│   ├── GALNT3, FUT8, ST6GAL1等关键酶分析
│   └── 酶-糖链调控网络可视化
├── 4. 多组学整合分析
│   ├── 糖组学 + 基因组学 + 蛋白质组学
│   ├── 联合生存分析
│   └── 通路富集分析（GSEA）
├── 5. 临床决策支持工具
│   ├── 列线图(Nomogram)在线生成
│   ├── 淋巴结转移风险预测
│   └── 患者分层与个性化报告
├── 6. 数据资源中心
│   ├── 公开数据集整合（GlycoPOST、TCGA等）
│   ├── 文献 curated 数据
│   └── Meta分析模块
└── 7. 知识库与文献
    ├── 糖基化-甲状腺癌关联知识图谱
    └── 自动文献更新
```

### 5.3 技术栈建议

| 层级 | 推荐技术 | 说明 |
|---|---|---|
| 前端 | React/Vue + D3.js/ECharts | 交互式可视化、糖链结构渲染 |
| 后端 | Python Flask/FastAPI | 生物信息学分析流程集成 |
| 数据库 | PostgreSQL + MongoDB | 结构化数据 + 糖谱非结构化数据 |
| 分析引擎 | R (ggplot2, pROC, caret) / Python (scikit-learn, pandas) | 统计分析和机器学习 |
| 糖组学工具 | GlycoWorkbench, SNFG符号 | 糖链结构表示 |
| 部署 | Docker + AWS/阿里云 | 容器化、可扩展 |

### 5.4 开发路线图

| 阶段 | 时间 | 目标 | 里程碑 |
|---|---|---|---|
| Phase 1 | 1-3个月 | 数据整合与基础数据库 | 整合已发表血清糖组学数据，建立基础数据库 |
| Phase 2 | 3-6个月 | 核心分析工具 | IgG N-糖组诊断模型在线验证工具上线 |
| Phase 3 | 6-12个月 | 扩展模块 | 组织糖组学、糖蛋白特异性分析、多组学整合 |
| Phase 4 | 12-18个月 | 临床转化 | 列线图生成、临床决策支持、外部验证 |

### 5.5 挑战与对策

| 挑战 | 对策 |
|---|---|
| 数据标准化问题（不同质谱平台、处理方法差异） | 建立统一数据预处理流程，支持多种格式导入 |
| 样本量限制（部分研究仅13-15例/组） | 整合多个研究公开数据，构建Meta分析模块 |
| 糖链结构注释复杂性 | 集成GlyTouCan和SNFG标准，自动化注释 |
| 临床验证需求 | 与医院合作，前瞻性队列验证 |

---

## 六、发表策略与影响力评估

### 6.1 目标期刊

| 期刊 | 类型 | 影响因子区间 | 适合阶段 |
|---|---|---|---|
| **Nucleic Acids Research** (Database Issue) | 数据库专刊 | ~14-15 | 平台v1.0完成后 |
| **Bioinformatics** | 方法学 | ~4-5 | 分析工具开发 |
| **Glycobiology** | 糖生物学专业 | ~3-4 | 糖组学特异性分析 |
| **Journal of Proteome Research** | 蛋白质组学 | ~4-5 | 糖蛋白组学整合 |
| **Frontiers in Oncology/Endocrinology** | 开放获取 | ~3-5 | 快速发表初步成果 |

### 6.2 预期影响力

- **科学价值**: 填补甲状腺癌糖基化领域专用平台空白
- **临床价值**: 为甲状腺癌诊断、预后评估提供无创糖基化标志物
- **方法学价值**: 建立糖组学-临床数据整合分析的标准化流程
- **社区价值**: 成为甲状腺癌糖基化研究的 central hub

---

## 附录：数据来源与脚本信息

### 自动化收集脚本

| 脚本 | 功能 | 输出 |
|---|---|---|
| `01_collect_literature.py` | 文献收集与整理 | `input/literature/` |
| `02_collect_databases.py` | 数据库信息收集 | `input/databases/` |
| `03_collect_datasets.py` | 数据集信息收集 | `input/datasets/` |
| `04_generate_report.py` | 综合报告生成 | `output/` |

### 输入数据结构
```
input/
├── literature/
│   ├── thyroid_cancer_glycosylation_literature.json
│   ├── thyroid_cancer_glycosylation_literature.csv
│   └── thyroid_cancer_glycosylation_literature.md
├── databases/
│   ├── glycosylation_databases.json
│   ├── glycosylation_databases.csv
│   └── glycosylation_databases.md
└── datasets/
    ├── thyroid_cancer_glycosylation_datasets.json
    ├── thyroid_cancer_glycosylation_datasets.csv
    └── thyroid_cancer_glycosylation_datasets.md
```

---

*本报告由自动化脚本生成，数据截至收集时间点。建议定期更新文献和数据集信息。*