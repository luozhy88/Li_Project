# Phase 1 交付报告：甲状腺癌糖基化基础数据库

> 生成时间: 2026-05-01 14:50
> 阶段目标: 整合已发表血清糖组学数据，建立基础数据库

---

## 一、Phase 1 执行摘要

### 1.1 交付成果

| 交付物 | 状态 | 说明 |
|---|---|---|
| SQLite基础数据库 | ✅ 完成 | `output/thyroid_glyco_db.sqlite` |
| 文献元数据整合 | ✅ 完成 | 11项研究，2018-2025 |
| 糖链结构库 | ✅ 完成 | 24种核心N-糖链结构 |
| 临床分组数据 | ✅ 完成 | 17个分组，涵盖HC/BTN/TC/复发等 |
| 糖链变化趋势 | ✅ 完成 | 21条丰度趋势记录（上调/下调） |
| 生物标志物库 | ✅ 完成 | 7个验证/候选标志物 |
| 糖基转移酶网络 | ✅ 完成 | 8个酶 + 12条酶-糖链关联 |
| 查询与浏览工具 | ✅ 完成 | `07_query_database.py` 交互式查询 |
| 数据导出功能 | ✅ 完成 | 支持CSV/JSON批量导出 |

### 1.2 数据库统计

| 数据类别 | 数量 |
|---|---|
| 研究/文献 | 11 项 |
| 糖链结构 | 24 种 |
| 临床分组 | 17 个 |
| 样本记录 | 17 条 |
| 糖链丰度趋势 | 21 条 |
| 生物标志物 | 7 个 |
| 糖基转移酶 | 8 个 |
| 酶-糖链关联 | 12 条 |

### 1.3 覆盖的癌种类型

- Differentiated thyroid cancer (DTC)
- Papillary thyroid carcinoma (PTC)
- Papillary thyroid microcarcinoma (PTMC)
- All types
- Medullary thyroid carcinoma (MTC)

### 1.4 研究方法分布

- MALDI-TOF MS
- HCD-MS/MS, MALDI-imaging
- 综述
- LC-MS/MS glycoproteomics
- LC-ESI/MS

---

## 二、数据库架构设计

采用 **SQLite** 关系型数据库，8个核心表：

```
thyroid_glyco_db.sqlite
├── studies                # 研究/文献元数据
├── glycan_structures      # 糖链结构信息 (SNFG, 分子量, 类型)
├── clinical_groups        # 临床分组 (HC, BTN, TC, 复发等)
├── samples                # 样本信息 (连接studies和groups)
├── glycan_abundance       # 糖链丰度/变化趋势 (核心数据)
├── biomarkers             # 已验证/候选生物标志物
├── glycosyltransferases   # 糖基转移酶信息
└── enzyme_glycan_links    # 酶-糖链调控网络
```

### ER关系图

```
studies (1) ────< clinical_groups (1) ────< samples (1) ────< glycan_abundance (N)
                │                                         ↑
                └── biomarkers (N)                          │
                                                            │
glycan_structures (1) <────── enzyme_glycan_links (N) ──────┘
                ↑
glycosyltransferases (1)
```

---

## 三、整合数据详情

### 3.1 研究覆盖度

| 研究ID | 文献标题 | 糖链趋势数据量 |
|---|---|---|
| S001 | Serum glycomic profile as a predictive biomarker o... | 6 |
| S002 | Diagnostic Potential of Plasma IgG N-glycans in Di... | 8 |
| S003 | Nomograms Based on Serum N-glycome for Diagnosis o... | 4 |
| S004 | Distinguishing Benign and Malignant Thyroid Nodule... | 1 |
| S005 | Thyroid Carcinoma Glycoproteins Express Altered N-... | 2 |
| S006 | Sialylation in Thyroid Carcinoma: An Overview of M... | 0 |
| S007 | N-glycan profiling of papillary thyroid carcinoma ... | 0 |
| S008 | Comparative glycoproteomic profiling of human body... | 0 |
| S009 | Glycosylation-Based Serum Biomarkers for Cancer Di... | 0 |
| S010 | Serum linkage-specific sialylation changes are pot... | 0 |
| S011 | Characterization of Human Medullary Thyroid Carcin... | 0 |

### 3.2 核心糖链结构库

已录入24种代表性N-糖链，覆盖以下类型：
- **Biantennary complex**: G0F, G1F, G2F, G0, G1, G2 及其唾液酸化变体
- **Bisecting type**: H3N5F1, H4N5F1, H5N5F1 (与IgG BN标志物相关)
- **High-mannose**: Man5-Man9 (TM标志物)
- **Multi-antennary**: CA3 (triantennary), CA4 (tetraantennary)
- **Special modifications**: 3-O-sulfated galactose, difucosylated
- **Hybrid type**: Hybrid structures

### 3.3 生物标志物清单

| 标志物 | 类型 | 样本 | AUC | 验证状态 |
|---|---|---|---|---|
| BN (bisecting type neutral N-g | 诊断标志物 | 血浆 IgG | 0.920 | 发现队列验证 |
| Glyco-panel | 诊断标志物 | 血浆 IgG | 0.917 | 发现队列+验证队列 |
| G0F:G1F ratio | 预后标志物 | 血清 | 0.820 | 发现队列+小样本验证 |
| TM + CA1 + CA4 + A2Fa | 诊断标志物 | 血清 | N/A | 单中心验证 |
| CA4 + A2F0S0G | 预后标志物 | 血清 | N/A | 单中心验证 |
| Plasma N-glycome panel | 诊断/预后标志物 | 血浆 | N/A | 发现阶段 |
| Linkage-specific sialylation | 预后标志物 | 血清 | N/A | 待验证 |

### 3.4 糖基转移酶网络

已整合8个关键糖基转移酶及其调控的糖链：
- **B4GALT1**: Beta-1,4-galactosyltransferase 1 (Beta-1,4 galactosylation)
- **FUT8**: Alpha-1,6-fucosyltransferase (Core fucosylation)
- **GAL3ST3**: Galactose-3-O-sulfotransferase 3 (3-O-sulfation)
- **GALNT3**: Polypeptide N-acetylgalactosaminyltransferase 3 (O-glycosylation initiation)
- **MAN1A2**: Mannosidase alpha class 1A member 2 (Mannose trimming)
- **MGAT5**: Alpha-1,6-mannosylglycoprotein 6-beta-N-acetylglucosaminyltransferase (Beta-1,6 branching)
- **ST3GAL**: ST3 beta-galactoside alpha-2,3-sialyltransferase (Alpha-2,3 sialylation)
- **ST6GAL1**: Beta-galactoside alpha-2,6-sialyltransferase 1 (Alpha-2,6 sialylation)

---

## 四、数据库使用指南

### 4.1 交互式查询

```bash
python3 07_query_database.py
```

支持以下查询模式：
1. **所有研究概览** - 快速浏览整合的文献
2. **特定研究糖链变化** - 输入研究ID查看该研究的糖链趋势
3. **生物标志物汇总** - 按AUC排序查看所有标志物
4. **糖链结构搜索** - 按组成或名称搜索糖链
5. **酶-糖链网络** - 查看糖基转移酶与糖链的调控关系
6. **跨研究比较** - 对比同一糖链在不同研究中的变化趋势

### 4.2 批量数据导出

```bash
python3 07_query_database.py --export
```

自动导出以下文件到 `output/queries/`：
- `all_studies.csv/json` - 所有研究
- `all_glycans.csv/json` - 所有糖链结构
- `all_biomarkers.csv/json` - 所有生物标志物
- `glycan_abundance_trends.csv/json` - 糖链丰度趋势
- `enzyme_glycan_network.csv/json` - 酶-糖链网络

### 4.3 直接SQL查询

```python
import sqlite3
conn = sqlite3.connect('output/thyroid_glyco_db.sqlite')
cursor = conn.cursor()
cursor.execute('SELECT * FROM biomarkers WHERE performance_auc > 0.8')
results = cursor.fetchall()
```

---

## 五、Phase 2 开发建议

基于Phase 1的基础数据库，建议Phase 2重点开发以下内容：

### 5.1 优先级1：IgG N-糖组诊断模型验证工具
- **目标**: 基于S002 (Zhang 2021) 的数据，开发在线BN特征计算和ROC分析工具
- **技术**: Python (scikit-learn) + 前端可视化
- **输入**: 用户上传的IgG N-糖组MALDI-TOF数据
- **输出**: 分类概率、ROC曲线、AUC报告

### 5.2 优先级2：数据量扩展
- 联系文献作者获取更完整的原始数据（丰度均值、标准差）
- 整合GlycoPOST公开数据集 (GPST000197, GPST000495)
- 导入TCGA糖基转移酶表达数据（RNAseq）

### 5.3 优先级3：Web界面原型
- 开发Flask/FastAPI后端API
- 实现基础的前端糖谱浏览页面
- 支持研究筛选、糖链搜索、结果可视化

---

## 附录：文件清单

### 脚本文件
| 脚本 | 功能 |
|---|---|
| `05_create_database.py` | 创建SQLite数据库架构 |
| `06_import_glycomics_data.py` | 导入文献糖组学数据 |
| `07_query_database.py` | 交互式查询与批量导出 |
| `08_generate_phase1_report.py` | 生成本报告 |

### 输出文件
| 文件 | 说明 |
|---|---|
| `output/thyroid_glyco_db.sqlite` | 核心数据库 |
| `output/phase1_report.md` | Phase 1交付报告 |
| `output/queries/*.csv` | 查询结果CSV |
| `output/queries/*.json` | 查询结果JSON |

---

*本报告由自动化脚本生成*