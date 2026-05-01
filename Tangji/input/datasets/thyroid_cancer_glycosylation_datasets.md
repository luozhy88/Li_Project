# 甲状腺癌糖基化相关可获取数据集

> 收集时间: 2026-05-01 14:32
> 数据集总数: 7 个
> 公开质谱数据: 2 个 (GlycoPOST)
> 文献衍生数据: 4 个
> TCGA基因组数据: 1 个

---

## DS001. GlycoPOST GPST000197

- **仓库**: GlycoPOST
- **URL**: https://glycopost.glycosmos.org/entry/GPST000197
- **登录号**: GPST000197
- **描述**: 人甲状腺髓样癌糖鞘脂的LC-ESI/MS原始数据
- **方法**: LC-ESI/MS (LTQ linear quadrupole ion trap)
- **数据格式**: 原始质谱文件 (.raw)
- **可用性**: 公开
- **MIRAGE合规**: 是
- **引用**: Characterization of Human Medullary Thyroid Carcinoma Glycosphingolipids (Cancers 2021)

**样本信息**:
- 类型: 组织
- 来源: Uppsala University Hospital, Sweden
- 疾病: Medullary thyroid carcinoma (MTC)
- 规模: 有限数量

---

## DS002. GlycoPOST GPST000495

- **仓库**: GlycoPOST
- **URL**: https://glycopost.glycosmos.org/entry/GPST000495
- **登录号**: GPST000495
- **描述**: 甲状腺癌糖蛋白异常N-聚糖（含3-O-硫酸化半乳糖残基）的质谱数据
- **方法**: HCD-MS/MS, MALDI-imaging
- **数据格式**: 原始质谱文件
- **可用性**: 公开
- **MIRAGE合规**: 是
- **引用**: Thyroid Carcinoma Glycoproteins Express Altered N-Glycans with 3-O-Sulfated Galactose Residues

**样本信息**:
- 类型: 组织
- 来源: Beth Israel Deaconess Medical Center
- 疾病: Papillary thyroid carcinoma (PTC)
- 规模: 配对PTC组织和正常组织

---

## DS003. TCGA-THCA RNAseq

- **仓库**: GDC Data Portal (TCGA)
- **URL**: https://portal.gdc.cancer.gov/
- **登录号**: TCGA-THCA
- **描述**: 甲状腺癌的mRNA转录组数据，可用于糖基转移酶表达分析
- **方法**: RNA sequencing
- **数据格式**: HTSeq-FPKM, counts
- **可用性**: 公开
- **MIRAGE合规**: N/A
- **引用**: TCGA Research Network
- **备注**: 虽非直接糖组学数据，但可用于糖基转移酶(GALNT3, FUT8, ST6GAL1等)的表达关联分析

**样本信息**:
- 类型: 组织
- 来源: 多中心收集
- 疾病: Thyroid carcinoma (THCA)
- 规模: 509例肿瘤组织，58例正常组织，共567例

---

## DS004. 血浆IgG N-糖组诊断数据集 (PUMCH)

- **仓库**: 文献原始数据
- **URL**: https://pmc.ncbi.nlm.nih.gov/articles/PMC8406750/
- **登录号**: N/A
- **描述**: 北京协和医院收集的血浆IgG N-糖组数据，用于区分甲状腺癌与良性结节
- **方法**: MALDI-TOF MS
- **数据格式**: 文献中提供的相对定量数据
- **可用性**: 需联系作者
- **MIRAGE合规**: 未知
- **引用**: Zhang ZJ et al. Front Oncol. 2021
- **备注**: IgG N-糖组诊断模型最成熟的标志物体系

**样本信息**:
- 类型: 血浆
- 来源: Peking Union Medical College Hospital
- 疾病: Papillary thyroid carcinoma (PTC), Benign thyroid nodules (BTN)
- 规模: 发现队列：TC(n=25) vs HC(n=25)；验证队列：TC(n=47), BTN(n=18), HC(n=44)

---

## DS005. 血清N-糖组复发预测数据集

- **仓库**: 文献原始数据
- **URL**: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10067050/
- **登录号**: N/A
- **描述**: 复发性分化型甲状腺癌患者血清N-糖组MALDI-TOF MS数据
- **方法**: MALDI-TOF MS (permethylation labeling)
- **数据格式**: 文献中提供的35种N-聚糖组成和相对丰度
- **可用性**: 需联系作者
- **MIRAGE合规**: 未知
- **引用**: Kudelka MR et al. Cancer Med. 2023
- **备注**: G0F:G1F比值作为复发预测标志物

**样本信息**:
- 类型: 血清
- 来源: 
- 疾病: Recurrent differentiated thyroid cancer
- 规模: 13例复发患者 vs 15例健康对照

---

## DS006. PTMC血清N-糖组数据集

- **仓库**: 文献原始数据
- **URL**: https://pmc.ncbi.nlm.nih.gov/articles/PMC9497917/
- **登录号**: N/A
- **描述**: 乳头状甲状腺微小癌(PTMC)血清N-糖组数据，含LNM预测
- **方法**: MALDI-TOF MS (Bruker rapifleXtreme)
- **数据格式**: 96种直接检测N-糖特征和91种衍生特征
- **可用性**: 需联系作者
- **MIRAGE合规**: 未知
- **引用**: Nomograms Based on Serum N-glycome for PTMC (Front Oncol. 2022)
- **备注**: 包含诊断和LNM预测的列线图模型

**样本信息**:
- 类型: 血清
- 来源: 
- 疾病: Papillary thyroid microcarcinoma (PTMC) with/without LNM
- 规模: 多组比较：PTMC+LNM, PTMC-LNM, HC

---

## DS007. 血浆N-糖组BTN/TC/LNM数据集

- **仓库**: 文献原始数据
- **URL**: https://pmc.ncbi.nlm.nih.gov/articles/PMC8267918/
- **登录号**: N/A
- **描述**: 血浆蛋白N-糖组数据，用于区分BTN和TC，并预测淋巴结转移
- **方法**: MALDI-TOF MS (含连接特异性唾液酸化信息)
- **数据格式**: 文献中提供
- **可用性**: 需联系作者
- **MIRAGE合规**: 未知
- **引用**: Front Endocrinol. 2021
- **备注**: 包含连接特异性唾液酸化分析

**样本信息**:
- 类型: 血浆
- 来源: 
- 疾病: Thyroid cancer (TC), Benign thyroid nodules (BTN)
- 规模: TC, BTN, HC三组

---
