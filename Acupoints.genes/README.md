# 针灸穴位-基因关联数据库 (Acupoint-Gene DB)

> 面向 SCI 发表的系统性文献数据库与分析平台

## 项目概述

本项目旨在构建一个系统性的**针灸穴位-基因关联数据库**，通过自动化采集 PubMed 文献，提取穴位-基因-疾病关联信息，为针灸分子机制研究提供数据支撑，最终服务于 SCI 论文发表。

## 项目结构

```
Acupoints.genes/
├── input/                  # 原始输入数据
│   ├── papers/            # 下载的PDF全文
│   └── xml/               # PubMed原始XML
├── output/                # 整理输出数据
│   ├── acupoint_gene.db   # SQLite主数据库
│   ├── *.csv              # 导出的表格
│   └── summary.json       # 汇总统计
├── src/                   # 源代码
│   ├── config.py          # 项目配置与检索策略
│   ├── fetcher/
│   │   └── pubmed_fetcher.py   # PubMed文献下载
│   ├── parser/
│   │   └── entity_extractor.py # 实体提取(穴位/基因/疾病)
│   ├── database/
│   │   └── models.py      # 数据库模型
│   ├── analyzer/
│   │   └── cli_analyzer.py     # CLI分析工具
│   └── web/
│       └── app.py         # Web分析平台
├── run_pipeline.py        # 主运行脚本
└── logs/                  # 运行日志
```

## 快速开始

### 1. 初始化数据库

```bash
python3 run_pipeline.py --init
```

### 2. 采集文献

```bash
# 采集300篇（默认）
python3 run_pipeline.py --fetch

# 或指定数量
python3 run_pipeline.py --fetch --target 500
```

### 3. 实体提取

```bash
python3 -m src.parser.entity_extractor
```

### 4. 启动分析平台

```bash
# CLI分析
python3 -m src.analyzer.cli_analyzer --all

# Web平台
python3 -m src.web.app --port 8080
# 访问 http://localhost:8080
```

### 5. 一键执行全部

```bash
python3 run_pipeline.py --all --target 300
python3 -m src.parser.entity_extractor
python3 -m src.analyzer.cli_analyzer --all
```

## 数据库Schema

### articles（文献主表）
- `pmid`, `pmcid`, `doi`, `title`, `abstract`
- `journal`, `year`, `authors`, `keywords`, `mesh_terms`
- `has_fulltext`, `pdf_path`, `xml_path`
- `search_strategy`

### genes（基因表）
- `gene_symbol`, `gene_name`, `entrez_id`

### acupoints（穴位表）
- `name_en`, `name_cn`, `code`, `meridian`

### article_gene / article_acupoint / article_disease（关联表）

## 检索策略

覆盖20+种检索策略，包括：
- 核心策略：`acupuncture + gene/expression`
- 穴位策略：太冲(LR3)、足三里(ST36)、百会(GV20)等
- 疾病策略：抑郁、高血压、帕金森、疼痛、炎症等
- 技术策略：microarray、RNA-seq、proteomics

## Web平台功能

- 📊 **数据概览**：文献、基因、穴位、疾病统计
- 📄 **文献列表**：支持搜索、分页
- 🔍 **文献详情**：摘要、关联实体、外部链接
- 🕸️ **共现网络**：穴位-基因关系可视化（ECharts）

## 输出文件

| 文件 | 说明 |
|------|------|
| `output/acupoint_gene.db` | SQLite数据库 |
| `output/articles.csv` | 文献清单 |
| `output/genes.csv` | 基因统计 |
| `output/acupoints.csv` | 穴位统计 |
| `output/summary.json` | 完整汇总JSON |
| `input/papers/*.pdf` | 开放获取全文 |
| `input/xml/*.xml` | PubMed原始XML |

## 依赖

- Python 3.8+
- jinja2（Web模板）
- 标准库：sqlite3, urllib, http.server, xml.etree

## 注意事项

1. **NCBI API限速**：无API Key约3 req/sec，建议申请NCBI API Key
2. **全文下载**：仅下载PMC开放获取（OA）文章
3. **实体识别**：基于词典匹配，后续可接入NLP模型提升精度

## 发表建议

基于本数据库可发表的SCI方向：
- 数据库论文（Database/Scientific Data）
- 系统综述/Meta分析（网络药理学框架）
- 特定穴位-基因关联的深入分析
- 多组学整合分析

## License

仅供学术研究使用。
