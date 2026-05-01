#!/usr/bin/env python3
"""
04_generate_report.py
生成甲状腺癌糖基化交互平台可行性分析综合报告，输出到 output/
"""

import json
import os
from datetime import datetime

INPUT_LIT = "input/literature/thyroid_cancer_glycosylation_literature.json"
INPUT_DB = "input/databases/glycosylation_databases.json"
INPUT_DS = "input/datasets/thyroid_cancer_glycosylation_datasets.json"
OUTPUT_DIR = "output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

def generate_markdown_report():
    lit = load_json(INPUT_LIT)
    db = load_json(INPUT_DB)
    ds = load_json(INPUT_DS)

    lit_list = lit["literature"]
    db_list = db["databases"]
    ds_list = ds["datasets"]

    # 统计
    categories = {}
    for item in lit_list:
        cat = item["category"]
        categories[cat] = categories.get(cat, 0) + 1

    methods = {}
    for item in lit_list:
        m = item["method"]
        methods[m] = methods.get(m, 0) + 1

    sample_types = {}
    for item in lit_list:
        s = item["sample_type"]
        sample_types[s] = sample_types.get(s, 0) + 1

    cancer_types = {}
    for item in lit_list:
        c = item["thyroid_cancer_type"]
        cancer_types[c] = cancer_types.get(c, 0) + 1

    lines = []
    lines.append("# 甲状腺癌糖基化交互分析平台 — 可行性分析综合报告")
    lines.append("")
    lines.append(f"> 报告生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    lines.append("> 基于自动化文献、数据库和数据集收集脚本生成")
    lines.append("")
    lines.append("---")
    lines.append("")

    # 摘要
    lines.append("## 一、执行摘要")
    lines.append("")
    lines.append("### 1.1 核心结论")
    lines.append("")
    lines.append("| 评估维度 | 评分 | 说明 |")
    lines.append("|---|---|---|")
    lines.append("| 科学意义 | ★★★★★ | 糖基化是甲状腺癌诊断、预后、治疗靶点的关键方向 |")
    lines.append("| 数据基础 | ★★★★☆ | 已有足够的发表数据可供整合和验证 |")
    lines.append("| 技术成熟度 | ★★★★★ | 糖组学质谱技术标准化程度高，MALDI-TOF MS等方法成熟 |")
    lines.append("| 临床需求 | ★★★★★ | 现有标志物（如Tg）存在明显局限性，亟需补充 |")
    lines.append("| 创新空间 | ★★★★★ | **尚无同类专用平台，领域空白** |")
    lines.append("| **总体可行性** | **高** | **建议优先启动** |")
    lines.append("")
    lines.append("### 1.2 关键发现")
    lines.append("")
    lines.append("- **文献证据充分**: 收集到14篇高质量研究（2018-2025），涵盖诊断、预后、机制三大方向")
    lines.append("- **数据库空白**: 10个主流糖组学/癌症数据库中，**无一个专门针对甲状腺癌糖基化**的交互分析平台")
    lines.append("- **数据可获取**: 2个公开GlycoPOST质谱数据集 + 4个文献数据集 + TCGA转录组数据可用于整合")
    lines.append("- **标志物成熟**: IgG N-糖组（BN特征，AUC=0.920）和血清G0F:G1F比值（AUC=0.82）已具备临床转化基础")
    lines.append("")
    lines.append("---")
    lines.append("")

    # 第二部分：文献分析
    lines.append("## 二、文献证据分析")
    lines.append("")
    lines.append(f"### 2.1 文献概览（共 {len(lit_list)} 篇）")
    lines.append("")
    lines.append("#### 按研究方向分类")
    lines.append("")
    lines.append("| 分类 | 数量 | 代表性文献 |")
    lines.append("|---|---|---|")
    for cat, cnt in sorted(categories.items(), key=lambda x: -x[1]):
        reps = [x["id"] for x in lit_list if x["category"] == cat][:2]
        lines.append(f"| {cat} | {cnt} | {', '.join(reps)} |")
    lines.append("")

    lines.append("#### 按样本类型分类")
    lines.append("")
    lines.append("| 样本类型 | 数量 |")
    lines.append("|---|---|")
    for st, cnt in sorted(sample_types.items(), key=lambda x: -x[1]):
        lines.append(f"| {st} | {cnt} |")
    lines.append("")

    lines.append("#### 按癌种类型分类")
    lines.append("")
    lines.append("| 癌种 | 数量 |")
    lines.append("|---|---|")
    for ct, cnt in sorted(cancer_types.items(), key=lambda x: -x[1]):
        lines.append(f"| {ct} | {cnt} |")
    lines.append("")

    lines.append("### 2.2 关键生物标志物汇总")
    lines.append("")
    lines.append("| 标志物 | 样本类型 | 性能 | 应用方向 | 文献 |")
    lines.append("|---|---|---|---|---|")
    for item in lit_list:
        if item.get("biomarker") and item["biomarker"] != "N/A":
            perf = ""
            for finding in item["key_findings"]:
                if "AUC" in finding:
                    perf = finding
                    break
            lines.append(f"| {item['biomarker']} | {item['sample_type']} | {perf or '见文献'} | {item['category']} | {item['id']} |")
    lines.append("")

    lines.append("### 2.3 核心文献详细列表")
    lines.append("")
    for item in lit_list:
        lines.append(f"#### {item['id']}: {item['title']}")
        lines.append(f"- **作者**: {item['authors']} | **年份**: {item['year']} | **期刊**: {item['journal']}")
        lines.append(f"- **PMID**: {item['pmid'] or 'N/A'} | **DOI**: {item['doi']}")
        lines.append(f"- **类型**: {item['category']} | **方法**: {item['method']} | **样本**: {item['sample_type']}")
        lines.append(f"- **标志物**: {item['biomarker']}")
        lines.append("- **主要发现**:")
        for finding in item["key_findings"]:
            lines.append(f"  - {finding}")
        lines.append("")
    lines.append("---")
    lines.append("")

    # 第三部分：数据库空白分析
    lines.append("## 三、现有数据库与平台空白分析")
    lines.append("")
    lines.append(f"### 3.1 已调查数据库（共 {len(db_list)} 个）")
    lines.append("")
    lines.append("| 数据库 | 类型 | 甲状腺癌特异性 | 糖基化数据 | URL |")
    lines.append("|---|---|---|---|---|")
    for item in db_list:
        tc_spec = "❌ 无" if "无特异性" in item["thyroid_cancer_data"] else "⚠️ 有限"
        glyco = "✅" if item["type"] in ["质谱原始数据存储库", "糖蛋白知识库", "色谱和质谱数据库", "国际糖链结构库", "糖蛋白数据库"] else "❌"
        lines.append(f"| {item['name']} | {item['type']} | {tc_spec} | {glyco} | [链接]({item['url']}) |")
    lines.append("")

    lines.append("### 3.2 领域空白总结")
    lines.append("")
    lines.append("**核心空白**: 目前尚无专门针对甲状腺癌糖基化的交互分析平台。")
    lines.append("")
    lines.append("现有资源的主要局限：")
    lines.append("1. **通用糖组学数据库**（如GlycoStore、GlycoPost）: 仅存储原始质谱数据，缺乏疾病特异性分析工具")
    lines.append("2. **癌症数据库**（如TCGA、THPA）: 以基因组/蛋白质组为主，糖基化数据整合度极低")
    lines.append("3. **糖蛋白数据库**（如UniCarbKB）: 聚焦结构信息，无临床决策支持功能")
    lines.append("4. **缺乏交互分析**: 无在线标志物筛选、机器学习验证、列线图生成等功能")
    lines.append("")
    lines.append("---")
    lines.append("")

    # 第四部分：数据集评估
    lines.append("## 四、可用数据集评估")
    lines.append("")
    lines.append(f"### 4.1 数据集清单（共 {len(ds_list)} 个）")
    lines.append("")
    lines.append("| 编号 | 名称 | 仓库 | 样本类型 | 疾病 | 可用性 |")
    lines.append("|---|---|---|---|---|---|")
    for item in ds_list:
        si = item.get("sample_info", {})
        lines.append(f"| {item['id']} | {item['name']} | {item['repository']} | {si.get('type','N/A')} | {si.get('disease','N/A')} | {item['availability']} |")
    lines.append("")

    lines.append("### 4.2 数据可用性评估")
    lines.append("")
    lines.append("| 数据类型 | 可用性 | 说明 |")
    lines.append("|---|---|---|")
    lines.append("| 血清/血浆糖组学数据 | ⭐⭐⭐⭐ 较丰富 | 多项研究已发表MALDI-TOF MS数据，涵盖PTC、PTMC、复发患者等 |")
    lines.append("| 组织糖组学数据 | ⭐⭐⭐ 中等 | 有FFPE组织和新鲜组织的N-糖谱数据（GlycoPOST） |")
    lines.append("| IgG糖基化数据 | ⭐⭐⭐⭐ 较丰富 | 作为系统性免疫标志物，已有大量队列数据 |")
    lines.append("| 糖基转移酶表达数据 | ⭐⭐ 有限 | 部分RNAseq和IHC数据可用（TCGA、THPA） |")
    lines.append("| 多组学整合数据 | ⭐⭐ 较少 | 蛋白质组学+糖蛋白组学联合分析的研究开始出现 |")
    lines.append("")
    lines.append("---")
    lines.append("")

    # 第五部分：平台设计建议
    lines.append("## 五、平台构建建议")
    lines.append("")
    lines.append("### 5.1 建议平台名称")
    lines.append("- **ThyGlycoPortal** (Thyroid Glycosylation Portal)")
    lines.append("- **TC-GlycoDB** (Thyroid Cancer Glycosylation Database)")
    lines.append("- **Glyco-Thyroid Atlas**")
    lines.append("")
    lines.append("### 5.2 核心功能模块")
    lines.append("")
    lines.append("```")
    lines.append("ThyGlycoPortal/")
    lines.append("├── 1. 糖谱浏览与比较模块")
    lines.append("│   ├── 正常 vs 癌组织/血清糖基化差异可视化")
    lines.append("│   ├── 不同亚型（PTC/FTC/MTC/ATC）糖谱比较")
    lines.append("│   └── 个体化糖谱上传与注释")
    lines.append("├── 2. 标志物筛选与验证工具")
    lines.append("│   ├── 基于机器学习的良恶性分类模型")
    lines.append("│   ├── 复发风险预测模型（整合G0F:G1F等）")
    lines.append("│   └── ROC曲线、AUC、敏感性/特异性计算")
    lines.append("├── 3. 糖基转移酶关联分析")
    lines.append("│   ├── 糖谱变化与酶表达关联（TCGA RNAseq）")
    lines.append("│   ├── GALNT3, FUT8, ST6GAL1等关键酶分析")
    lines.append("│   └── 酶-糖链调控网络可视化")
    lines.append("├── 4. 多组学整合分析")
    lines.append("│   ├── 糖组学 + 基因组学 + 蛋白质组学")
    lines.append("│   ├── 联合生存分析")
    lines.append("│   └── 通路富集分析（GSEA）")
    lines.append("├── 5. 临床决策支持工具")
    lines.append("│   ├── 列线图(Nomogram)在线生成")
    lines.append("│   ├── 淋巴结转移风险预测")
    lines.append("│   └── 患者分层与个性化报告")
    lines.append("├── 6. 数据资源中心")
    lines.append("│   ├── 公开数据集整合（GlycoPOST、TCGA等）")
    lines.append("│   ├── 文献 curated 数据")
    lines.append("│   └── Meta分析模块")
    lines.append("└── 7. 知识库与文献")
    lines.append("    ├── 糖基化-甲状腺癌关联知识图谱")
    lines.append("    └── 自动文献更新")
    lines.append("```")
    lines.append("")

    lines.append("### 5.3 技术栈建议")
    lines.append("")
    lines.append("| 层级 | 推荐技术 | 说明 |")
    lines.append("|---|---|---|")
    lines.append("| 前端 | React/Vue + D3.js/ECharts | 交互式可视化、糖链结构渲染 |")
    lines.append("| 后端 | Python Flask/FastAPI | 生物信息学分析流程集成 |")
    lines.append("| 数据库 | PostgreSQL + MongoDB | 结构化数据 + 糖谱非结构化数据 |")
    lines.append("| 分析引擎 | R (ggplot2, pROC, caret) / Python (scikit-learn, pandas) | 统计分析和机器学习 |")
    lines.append("| 糖组学工具 | GlycoWorkbench, SNFG符号 | 糖链结构表示 |")
    lines.append("| 部署 | Docker + AWS/阿里云 | 容器化、可扩展 |")
    lines.append("")

    lines.append("### 5.4 开发路线图")
    lines.append("")
    lines.append("| 阶段 | 时间 | 目标 | 里程碑 |")
    lines.append("|---|---|---|---|")
    lines.append("| Phase 1 | 1-3个月 | 数据整合与基础数据库 | 整合已发表血清糖组学数据，建立基础数据库 |")
    lines.append("| Phase 2 | 3-6个月 | 核心分析工具 | IgG N-糖组诊断模型在线验证工具上线 |")
    lines.append("| Phase 3 | 6-12个月 | 扩展模块 | 组织糖组学、糖蛋白特异性分析、多组学整合 |")
    lines.append("| Phase 4 | 12-18个月 | 临床转化 | 列线图生成、临床决策支持、外部验证 |")
    lines.append("")

    lines.append("### 5.5 挑战与对策")
    lines.append("")
    lines.append("| 挑战 | 对策 |")
    lines.append("|---|---|")
    lines.append("| 数据标准化问题（不同质谱平台、处理方法差异） | 建立统一数据预处理流程，支持多种格式导入 |")
    lines.append("| 样本量限制（部分研究仅13-15例/组） | 整合多个研究公开数据，构建Meta分析模块 |")
    lines.append("| 糖链结构注释复杂性 | 集成GlyTouCan和SNFG标准，自动化注释 |")
    lines.append("| 临床验证需求 | 与医院合作，前瞻性队列验证 |")
    lines.append("")
    lines.append("---")
    lines.append("")

    # 第六部分：发表策略
    lines.append("## 六、发表策略与影响力评估")
    lines.append("")
    lines.append("### 6.1 目标期刊")
    lines.append("")
    lines.append("| 期刊 | 类型 | 影响因子区间 | 适合阶段 |")
    lines.append("|---|---|---|---|")
    lines.append("| **Nucleic Acids Research** (Database Issue) | 数据库专刊 | ~14-15 | 平台v1.0完成后 |")
    lines.append("| **Bioinformatics** | 方法学 | ~4-5 | 分析工具开发 |")
    lines.append("| **Glycobiology** | 糖生物学专业 | ~3-4 | 糖组学特异性分析 |")
    lines.append("| **Journal of Proteome Research** | 蛋白质组学 | ~4-5 | 糖蛋白组学整合 |")
    lines.append("| **Frontiers in Oncology/Endocrinology** | 开放获取 | ~3-5 | 快速发表初步成果 |")
    lines.append("")

    lines.append("### 6.2 预期影响力")
    lines.append("")
    lines.append("- **科学价值**: 填补甲状腺癌糖基化领域专用平台空白")
    lines.append("- **临床价值**: 为甲状腺癌诊断、预后评估提供无创糖基化标志物")
    lines.append("- **方法学价值**: 建立糖组学-临床数据整合分析的标准化流程")
    lines.append("- **社区价值**: 成为甲状腺癌糖基化研究的 central hub")
    lines.append("")
    lines.append("---")
    lines.append("")

    # 附录
    lines.append("## 附录：数据来源与脚本信息")
    lines.append("")
    lines.append("### 自动化收集脚本")
    lines.append("")
    lines.append("| 脚本 | 功能 | 输出 |")
    lines.append("|---|---|---|")
    lines.append("| `01_collect_literature.py` | 文献收集与整理 | `input/literature/` |")
    lines.append("| `02_collect_databases.py` | 数据库信息收集 | `input/databases/` |")
    lines.append("| `03_collect_datasets.py` | 数据集信息收集 | `input/datasets/` |")
    lines.append("| `04_generate_report.py` | 综合报告生成 | `output/` |")
    lines.append("")
    lines.append("### 输入数据结构")
    lines.append("```")
    lines.append("input/")
    lines.append("├── literature/")
    lines.append("│   ├── thyroid_cancer_glycosylation_literature.json")
    lines.append("│   ├── thyroid_cancer_glycosylation_literature.csv")
    lines.append("│   └── thyroid_cancer_glycosylation_literature.md")
    lines.append("├── databases/")
    lines.append("│   ├── glycosylation_databases.json")
    lines.append("│   ├── glycosylation_databases.csv")
    lines.append("│   └── glycosylation_databases.md")
    lines.append("└── datasets/")
    lines.append("    ├── thyroid_cancer_glycosylation_datasets.json")
    lines.append("    ├── thyroid_cancer_glycosylation_datasets.csv")
    lines.append("    └── thyroid_cancer_glycosylation_datasets.md")
    lines.append("```")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("*本报告由自动化脚本生成，数据截至收集时间点。建议定期更新文献和数据集信息。*")

    path = os.path.join(OUTPUT_DIR, "feasibility_report.md")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"[OK] Markdown report saved: {path}")

def generate_json_summary():
    lit = load_json(INPUT_LIT)
    db = load_json(INPUT_DB)
    ds = load_json(INPUT_DS)

    summary = {
        "metadata": {
            "title": "甲状腺癌糖基化交互分析平台可行性分析报告",
            "generated": datetime.now().isoformat(),
            "version": "1.0"
        },
        "executive_summary": {
            "overall_feasibility": "高",
            "scientific_significance": "★★★★★",
            "data_foundation": "★★★★☆",
            "technical_maturity": "★★★★★",
            "clinical_need": "★★★★★",
            "innovation_space": "★★★★★",
            "key_finding": "目前尚无专门针对甲状腺癌糖基化的交互分析平台，领域存在重大空白"
        },
        "statistics": {
            "literature_count": len(lit["literature"]),
            "database_count": len(db["databases"]),
            "dataset_count": len(ds["datasets"]),
            "categories": list(set(x["category"] for x in lit["literature"])),
            "database_gap": "无甲状腺癌糖基化专用平台"
        },
        "recommended_platform_name": ["ThyGlycoPortal", "TC-GlycoDB", "Glyco-Thyroid Atlas"],
        "core_modules": [
            "糖谱浏览与比较模块",
            "标志物筛选与验证工具",
            "糖基转移酶关联分析",
            "多组学整合分析",
            "临床决策支持工具",
            "数据资源中心",
            "知识库与文献"
        ],
        "publication_targets": [
            {"journal": "Nucleic Acids Research (Database Issue)", "timing": "平台v1.0完成后"},
            {"journal": "Bioinformatics", "timing": "分析工具开发阶段"},
            {"journal": "Glycobiology", "timing": "糖组学特异性分析"}
        ]
    }

    path = os.path.join(OUTPUT_DIR, "feasibility_report.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    print(f"[OK] JSON summary saved: {path}")

if __name__ == "__main__":
    generate_markdown_report()
    generate_json_summary()
    print("\n综合报告生成完成。")
