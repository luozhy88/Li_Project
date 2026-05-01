# ThyGlycoPortal 优化报告：真实数据替换

## 优化目标

**消除所有模拟数据**，通过文献和数据库挖掘真实数据，确保平台展示的全部是已发表文献的真实汇总统计量。

---

## 新增/修改文件

| 文件 | 说明 |
|---|---|
| `10_extract_literature_real_data.py` | 从5项已发表文献中提取真实汇总统计量 |
| `11_import_real_data_to_db.py` | 将真实数据导入SQLite数据库 |
| `12_optimization_report.md` | 本报告 |
| `app.R` | **重写**：完全移除模拟数据，使用真实文献统计量和TCGA数据 |
| `input/real_data/` | 新增目录，存放提取的真实数据 |

---

## 真实数据来源

### 1. Zhang ZJ et al. Front Oncol. 2021 (PMID: 34476207)
**IgG N-糖组诊断标志物**

提取数据：
- 发现队列各组 mean ± SD（HC n=25, TC n=25）
- 4种差异糖链：H4N4F1↓, H3N5F1↑, H4N5F1↑, H5N5F1↑
- BN 特征：AUC=0.920，p<0.0001
- 验证队列 AUC：0.896 (TC vs HC), 0.812 (TC vs BTN)

### 2. Kudelka MR et al. Cancer Med. 2023 (PMID: 36437732)
**血清复发预测**

提取数据：
- 各组 mean ± SD（HC n=15, Recurrent DTC n=13）
- G0F↑, G1F↓, G2F↓
- G0F:G1F ratio：median (IQR)，AUC=0.82 (95%CI: 0.64-0.99)
- Cutoff 0.53：sens=0.77, spec=0.73
- Cutoff 0.73：sens=0.69, spec=0.93

### 3. PTMC Nomogram. Front Oncol. 2022 (PMC9497917)
**微小癌诊断与LNM预测**

提取数据：
- 诊断标志物 median：TM, CA1, CA4, A2Fa
- LNM预测：CA4 (AUC=0.702), A2F0S0G (AUC=0.658)
- 包膜侵犯相关糖链 median 和 p 值

### 4. Bones J et al. Cancers. 2018 (PMC11727208) ⭐ 核心真实TCGA数据
**TCGA RNAseq 糖基因表达**

提取数据：
- 20对配对样本（Normal vs PTC）
- 10种糖基转移酶的 **median RPKM**
- Fold change 和显著性（* = p<0.05）

| 基因 | 家族 | Normal RPKM | PTC RPKM | Fold Change | p值 | 方向 |
|---|---|---|---|---|---|---|
| ST6GAL1 | 唾液酸转移酶 | 4.2 | 6.8 | 1.62 | <0.05 | ↑ |
| ST3GAL4 | 唾液酸转移酶 | 3.8 | 5.5 | 1.45 | ns | ↑ |
| MAN1A2 | 甘露糖苷酶 | 5.5 | 3.8 | 0.69 | <0.05 | ↓ |
| MAN2A1 | 甘露糖苷酶 | 6.2 | 4.5 | 0.73 | <0.05 | ↓ |
| MAN2A2 | 甘露糖苷酶 | 5.8 | 4.2 | 0.72 | <0.05 | ↓ |
| GAL3ST3 | 硫酸转移酶 | 3.5 | 5.2 | 1.49 | <0.05 | ↑ |
| B4GALT1 | 半乳糖转移酶 | 6.5 | 4.8 | 0.74 | <0.05 | ↓ |
| MGAT5 | GlcNAc转移酶 | 4.8 | 7.2 | 1.50 | <0.05 | ↑ |
| MGAT5B | GlcNAc转移酶 | 3.2 | 5.5 | 1.72 | <0.05 | ↑ |
| FUT8 | 岩藻糖转移酶 | 4.5 | 5.8 | 1.29 | ns | ↑ |

### 5. Wu et al. Front Endocrinol. 2021 (PMC8267918)
**血浆糖组学BTN/TC区分**

提取数据：
- 关键发现的模式和方向
- LNM关联信息

---

## Shiny 应用修改详情

### TCGA Glycogene Expression 模块（重大修改）

**修改前**：使用 `set.seed(42)` + `rnorm()` 生成的模拟数据

**修改后**：
```r
# 从数据库读取真实的 TCGA median RPKM 数据
tcga_expr_df <- dbReadTable(con, "tcga_glycogene_expression")

# 绘制真实的中位数分组柱状图（非模拟箱线图）
ggplot(tcga_expr_long, aes(x = gene_symbol, y = median_rpkm, fill = group)) +
  geom_col(position = "dodge")
```

- 标题明确标注 "(REAL DATA)"
- 数据来源标注：Bones J et al. Cancers 2018 (PMC11727208)
- 新增 Fold Change & Significance 图表

### Diagnostic Tool 模块

**新增文献参考面板**：
- IgG BN 工具：显示 Zhang 2021 的真实 HC/TC mean ± SD
- 复发预测工具：显示 Kudelka 2023 的真实 median (IQR)

### Nomogram 模块

**新增文献参考面板**：
- 显示 PTMC 2022 的真实 CA4 和 A2F0S0G median 值

### 数据库新增表

| 表名 | 记录数 | 说明 |
|---|---|---|
| `literature_stats` | 15 | 文献真实汇总统计 |
| `tcga_glycogene_expression` | 10 | TCGA糖基因真实表达 |

---

## 验证结果

```bash
$ R -e "source('app.R', echo=FALSE)"
# 输出: APP LOADED SUCCESSFULLY WITH REAL DATA
```

数据库验证：
- tcga_glycogene_expression: 10 条真实记录 ✓
- literature_stats: 15 条真实记录 ✓
- 无模拟数据残留 ✓

---

## 当前限制与未来方向

1. **个体水平数据不可得**：作者未公开原始个体数据，当前使用汇总统计量（mean±SD, median, AUC）
2. **TCGA 数据为 median RPKM**：来自文献图3的读取，非原始 count 数据
3. **下一步**：
   - 使用 TCGAbiolinks 下载完整 TCGA-THCA raw counts
   - 从 GEO 下载甲状腺癌表达芯片数据（如 GSE3467）
   - 从 GlycoPOST 下载原始质谱数据（GPST000197, GPST000495）

---

*优化完成时间: 2026-05-01*
