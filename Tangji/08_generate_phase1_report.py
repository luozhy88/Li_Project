#!/usr/bin/env python3
"""
08_generate_phase1_report.py
Phase 1 交付报告生成：数据库建设、数据整合成果与下一步建议
"""

import sqlite3
import json
import os
from datetime import datetime

DB_PATH = "output/thyroid_glyco_db.sqlite"
OUTPUT_DIR = "output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def get_stats():
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    stats = {}

    for table in ["studies", "glycan_structures", "clinical_groups", "samples",
                  "glycan_abundance", "biomarkers", "glycosyltransferases", "enzyme_glycan_links"]:
        cursor.execute(f"SELECT COUNT(*) FROM {table}")
        stats[table] = cursor.fetchone()[0]

    cursor.execute("SELECT DISTINCT cancer_type FROM studies")
    stats["cancer_types"] = [r[0] for r in cursor.fetchall()]

    cursor.execute("SELECT DISTINCT method FROM studies")
    stats["methods"] = [r[0] for r in cursor.fetchall()]

    cursor.execute("SELECT MIN(year), MAX(year) FROM studies")
    stats["year_range"] = cursor.fetchone()

    cursor.execute("""
        SELECT s.study_id, s.title, COUNT(ga.abundance_id) as glycan_count
        FROM studies s
        LEFT JOIN samples sam ON s.study_id = sam.study_id
        LEFT JOIN glycan_abundance ga ON sam.sample_id = ga.sample_id
        GROUP BY s.study_id
    """)
    stats["study_glycan_coverage"] = cursor.fetchall()

    cursor.execute("""
        SELECT biomarker_type, COUNT(*) FROM biomarkers GROUP BY biomarker_type
    """)
    stats["biomarker_types"] = cursor.fetchall()

    conn.close()
    return stats

def generate_markdown_report():
    stats = get_stats()

    lines = []
    lines.append("# Phase 1 交付报告：甲状腺癌糖基化基础数据库")
    lines.append("")
    lines.append(f"> 生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    lines.append("> 阶段目标: 整合已发表血清糖组学数据，建立基础数据库")
    lines.append("")
    lines.append("---")
    lines.append("")

    # 1. 执行摘要
    lines.append("## 一、Phase 1 执行摘要")
    lines.append("")
    lines.append("### 1.1 交付成果")
    lines.append("")
    lines.append("| 交付物 | 状态 | 说明 |")
    lines.append("|---|---|---|")
    lines.append("| SQLite基础数据库 | ✅ 完成 | `output/thyroid_glyco_db.sqlite` |")
    lines.append("| 文献元数据整合 | ✅ 完成 | 11项研究，2018-2025 |")
    lines.append("| 糖链结构库 | ✅ 完成 | 24种核心N-糖链结构 |")
    lines.append("| 临床分组数据 | ✅ 完成 | 17个分组，涵盖HC/BTN/TC/复发等 |")
    lines.append("| 糖链变化趋势 | ✅ 完成 | 21条丰度趋势记录（上调/下调） |")
    lines.append("| 生物标志物库 | ✅ 完成 | 7个验证/候选标志物 |")
    lines.append("| 糖基转移酶网络 | ✅ 完成 | 8个酶 + 12条酶-糖链关联 |")
    lines.append("| 查询与浏览工具 | ✅ 完成 | `07_query_database.py` 交互式查询 |")
    lines.append("| 数据导出功能 | ✅ 完成 | 支持CSV/JSON批量导出 |")
    lines.append("")

    lines.append("### 1.2 数据库统计")
    lines.append("")
    lines.append("| 数据类别 | 数量 |")
    lines.append("|---|---|")
    lines.append(f"| 研究/文献 | {stats['studies']} 项 |")
    lines.append(f"| 糖链结构 | {stats['glycan_structures']} 种 |")
    lines.append(f"| 临床分组 | {stats['clinical_groups']} 个 |")
    lines.append(f"| 样本记录 | {stats['samples']} 条 |")
    lines.append(f"| 糖链丰度趋势 | {stats['glycan_abundance']} 条 |")
    lines.append(f"| 生物标志物 | {stats['biomarkers']} 个 |")
    lines.append(f"| 糖基转移酶 | {stats['glycosyltransferases']} 个 |")
    lines.append(f"| 酶-糖链关联 | {stats['enzyme_glycan_links']} 条 |")
    lines.append("")

    lines.append("### 1.3 覆盖的癌种类型")
    lines.append("")
    for ct in stats["cancer_types"]:
        lines.append(f"- {ct}")
    lines.append("")

    lines.append("### 1.4 研究方法分布")
    lines.append("")
    for m in stats["methods"]:
        lines.append(f"- {m}")
    lines.append("")
    lines.append("---")
    lines.append("")

    # 2. 数据库架构
    lines.append("## 二、数据库架构设计")
    lines.append("")
    lines.append("采用 **SQLite** 关系型数据库，8个核心表：")
    lines.append("")
    lines.append("```")
    lines.append("thyroid_glyco_db.sqlite")
    lines.append("├── studies                # 研究/文献元数据")
    lines.append("├── glycan_structures      # 糖链结构信息 (SNFG, 分子量, 类型)")
    lines.append("├── clinical_groups        # 临床分组 (HC, BTN, TC, 复发等)")
    lines.append("├── samples                # 样本信息 (连接studies和groups)")
    lines.append("├── glycan_abundance       # 糖链丰度/变化趋势 (核心数据)")
    lines.append("├── biomarkers             # 已验证/候选生物标志物")
    lines.append("├── glycosyltransferases   # 糖基转移酶信息")
    lines.append("└── enzyme_glycan_links    # 酶-糖链调控网络")
    lines.append("```")
    lines.append("")
    lines.append("### ER关系图")
    lines.append("")
    lines.append("```")
    lines.append("studies (1) ────< clinical_groups (1) ────< samples (1) ────< glycan_abundance (N)")
    lines.append("                │                                         ↑")
    lines.append("                └── biomarkers (N)                          │")
    lines.append("                                                            │")
    lines.append("glycan_structures (1) <────── enzyme_glycan_links (N) ──────┘")
    lines.append("                ↑")
    lines.append("glycosyltransferases (1)")
    lines.append("```")
    lines.append("")
    lines.append("---")
    lines.append("")

    # 3. 数据内容详情
    lines.append("## 三、整合数据详情")
    lines.append("")

    lines.append("### 3.1 研究覆盖度")
    lines.append("")
    lines.append("| 研究ID | 文献标题 | 糖链趋势数据量 |")
    lines.append("|---|---|---|")
    for study_id, title, count in stats["study_glycan_coverage"]:
        lines.append(f"| {study_id} | {title[:50]}... | {count} |")
    lines.append("")

    lines.append("### 3.2 核心糖链结构库")
    lines.append("")
    lines.append("已录入24种代表性N-糖链，覆盖以下类型：")
    lines.append("- **Biantennary complex**: G0F, G1F, G2F, G0, G1, G2 及其唾液酸化变体")
    lines.append("- **Bisecting type**: H3N5F1, H4N5F1, H5N5F1 (与IgG BN标志物相关)")
    lines.append("- **High-mannose**: Man5-Man9 (TM标志物)")
    lines.append("- **Multi-antennary**: CA3 (triantennary), CA4 (tetraantennary)")
    lines.append("- **Special modifications**: 3-O-sulfated galactose, difucosylated")
    lines.append("- **Hybrid type**: Hybrid structures")
    lines.append("")

    lines.append("### 3.3 生物标志物清单")
    lines.append("")
    lines.append("| 标志物 | 类型 | 样本 | AUC | 验证状态 |")
    lines.append("|---|---|---|---|---|")

    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()
    cursor.execute("""
        SELECT b.name, b.biomarker_type, b.sample_type, b.performance_auc, b.validation_status
        FROM biomarkers b ORDER BY b.performance_auc DESC
    """)
    for row in cursor.fetchall():
        auc = f"{row[3]:.3f}" if row[3] else "N/A"
        lines.append(f"| {row[0][:30]} | {row[1]} | {row[2]} | {auc} | {row[4]} |")
    lines.append("")

    lines.append("### 3.4 糖基转移酶网络")
    lines.append("")
    lines.append("已整合8个关键糖基转移酶及其调控的糖链：")
    cursor.execute("SELECT gene_symbol, enzyme_name, pathway FROM glycosyltransferases ORDER BY gene_symbol")
    for row in cursor.fetchall():
        lines.append(f"- **{row[0]}**: {row[1]} ({row[2]})")
    lines.append("")
    conn.close()

    lines.append("---")
    lines.append("")

    # 4. 使用指南
    lines.append("## 四、数据库使用指南")
    lines.append("")
    lines.append("### 4.1 交互式查询")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 07_query_database.py")
    lines.append("```")
    lines.append("")
    lines.append("支持以下查询模式：")
    lines.append("1. **所有研究概览** - 快速浏览整合的文献")
    lines.append("2. **特定研究糖链变化** - 输入研究ID查看该研究的糖链趋势")
    lines.append("3. **生物标志物汇总** - 按AUC排序查看所有标志物")
    lines.append("4. **糖链结构搜索** - 按组成或名称搜索糖链")
    lines.append("5. **酶-糖链网络** - 查看糖基转移酶与糖链的调控关系")
    lines.append("6. **跨研究比较** - 对比同一糖链在不同研究中的变化趋势")
    lines.append("")

    lines.append("### 4.2 批量数据导出")
    lines.append("")
    lines.append("```bash")
    lines.append("python3 07_query_database.py --export")
    lines.append("```")
    lines.append("")
    lines.append("自动导出以下文件到 `output/queries/`：")
    lines.append("- `all_studies.csv/json` - 所有研究")
    lines.append("- `all_glycans.csv/json` - 所有糖链结构")
    lines.append("- `all_biomarkers.csv/json` - 所有生物标志物")
    lines.append("- `glycan_abundance_trends.csv/json` - 糖链丰度趋势")
    lines.append("- `enzyme_glycan_network.csv/json` - 酶-糖链网络")
    lines.append("")

    lines.append("### 4.3 直接SQL查询")
    lines.append("")
    lines.append("```python")
    lines.append("import sqlite3")
    lines.append("conn = sqlite3.connect('output/thyroid_glyco_db.sqlite')")
    lines.append("cursor = conn.cursor()")
    lines.append("cursor.execute('SELECT * FROM biomarkers WHERE performance_auc > 0.8')")
    lines.append("results = cursor.fetchall()")
    lines.append("```")
    lines.append("")
    lines.append("---")
    lines.append("")

    # 5. Phase 2 建议
    lines.append("## 五、Phase 2 开发建议")
    lines.append("")
    lines.append("基于Phase 1的基础数据库，建议Phase 2重点开发以下内容：")
    lines.append("")
    lines.append("### 5.1 优先级1：IgG N-糖组诊断模型验证工具")
    lines.append("- **目标**: 基于S002 (Zhang 2021) 的数据，开发在线BN特征计算和ROC分析工具")
    lines.append("- **技术**: Python (scikit-learn) + 前端可视化")
    lines.append("- **输入**: 用户上传的IgG N-糖组MALDI-TOF数据")
    lines.append("- **输出**: 分类概率、ROC曲线、AUC报告")
    lines.append("")
    lines.append("### 5.2 优先级2：数据量扩展")
    lines.append("- 联系文献作者获取更完整的原始数据（丰度均值、标准差）")
    lines.append("- 整合GlycoPOST公开数据集 (GPST000197, GPST000495)")
    lines.append("- 导入TCGA糖基转移酶表达数据（RNAseq）")
    lines.append("")
    lines.append("### 5.3 优先级3：Web界面原型")
    lines.append("- 开发Flask/FastAPI后端API")
    lines.append("- 实现基础的前端糖谱浏览页面")
    lines.append("- 支持研究筛选、糖链搜索、结果可视化")
    lines.append("")
    lines.append("---")
    lines.append("")

    lines.append("## 附录：文件清单")
    lines.append("")
    lines.append("### 脚本文件")
    lines.append("| 脚本 | 功能 |")
    lines.append("|---|---|")
    lines.append("| `05_create_database.py` | 创建SQLite数据库架构 |")
    lines.append("| `06_import_glycomics_data.py` | 导入文献糖组学数据 |")
    lines.append("| `07_query_database.py` | 交互式查询与批量导出 |")
    lines.append("| `08_generate_phase1_report.py` | 生成本报告 |")
    lines.append("")
    lines.append("### 输出文件")
    lines.append("| 文件 | 说明 |")
    lines.append("|---|---|")
    lines.append("| `output/thyroid_glyco_db.sqlite` | 核心数据库 |")
    lines.append("| `output/phase1_report.md` | Phase 1交付报告 |")
    lines.append("| `output/queries/*.csv` | 查询结果CSV |")
    lines.append("| `output/queries/*.json` | 查询结果JSON |")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("*本报告由自动化脚本生成*")

    path = os.path.join(OUTPUT_DIR, "phase1_report.md")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"[OK] Phase 1报告已生成: {path}")

def generate_json_summary():
    stats = get_stats()
    summary = {
        "phase": "Phase 1",
        "title": "甲状腺癌糖基化基础数据库",
        "generated": datetime.now().isoformat(),
        "deliverables": {
            "database": "output/thyroid_glyco_db.sqlite",
            "report": "output/phase1_report.md",
            "query_tool": "07_query_database.py"
        },
        "statistics": {
            "studies": stats["studies"],
            "glycan_structures": stats["glycan_structures"],
            "clinical_groups": stats["clinical_groups"],
            "samples": stats["samples"],
            "glycan_abundance_records": stats["glycan_abundance"],
            "biomarkers": stats["biomarkers"],
            "glycosyltransferases": stats["glycosyltransferases"],
            "enzyme_glycan_links": stats["enzyme_glycan_links"]
        },
        "coverage": {
            "cancer_types": stats["cancer_types"],
            "methods": stats["methods"],
            "year_range": stats["year_range"]
        },
        "next_phase_priorities": [
            "IgG N-glycan diagnostic model validation tool",
            "Data volume expansion (contact authors, GlycoPOST)",
            "Web interface prototype (Flask + frontend)"
        ]
    }
    path = os.path.join(OUTPUT_DIR, "phase1_summary.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    print(f"[OK] Phase 1摘要JSON: {path}")

if __name__ == "__main__":
    print("=" * 60)
    print("生成 Phase 1 交付报告")
    print("=" * 60)
    generate_markdown_report()
    generate_json_summary()
    print("\nPhase 1 报告生成完成！")
