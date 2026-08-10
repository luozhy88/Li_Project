# 假定位体征（False Localizing Sign）文献计量学分析流程

本仓库用于对“假定位体征 / false localizing sign / false localising sign / pseudo-localization sign”相关文献进行可复现的文献计量学分析，并自动生成中文 Word 报告（图表英文展示）。

## 项目结构

```text
.
├── config.yaml                  # 配置文件：路径、检索式、分析参数等
├── requirements.txt             # Python 依赖
├── 0_run_all.py                 # 一键运行完整流程
├── 1_parse_and_merge.py         # 解析 PubMed / WoS、去重、导出 CiteSpace 格式
├── 2_descriptive_analysis.py    # 描述性统计与图表（年度/期刊/国家/作者）
├── 3_keyword_analysis.py        # 关键词词频、共现网络、词云
├── 4_generate_report.py         # 生成中文 Word 报告
├── input/                       # 输入数据（用户提供）
│   └── false/                   # PubMed 检索结果示例
└── output/                      # 所有输出文件
    ├── 01_merged_records.csv
    ├── 01_citespace_input.txt
    ├── figures/                 # 英文图表（PNG）
    ├── tables/                  # 统计表格（CSV）
    ├── report/                  # 中文 Word 报告
    └── logs/                    # 运行日志
```

## 快速开始

### 1. 创建虚拟环境并安装依赖

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### 2. 配置数据与参数

编辑 `config.yaml`：

- `input.pubmed_medline`：PubMed MEDLINE 格式文件路径。
- `input.wos_file`：Web of Science 纯文本/CiteSpace 格式文件路径（可选；如未提供则仅分析 PubMed）。
- `output.root`：输出根目录（默认 `output`）。
- `search`：检索词与检索式（用于报告与方法部分）。
- `analysis`：图表与分析参数（Top N、共现阈值等）。

### 3. 运行分析

一键运行：

```bash
python 0_run_all.py
```

或分步运行：

```bash
python 1_parse_and_merge.py
python 2_descriptive_analysis.py
python 3_keyword_analysis.py
python 4_generate_report.py
```

### 4. 查看日志

使用一键脚本运行时，所有步骤的完整终端输出会自动保存到本地日志，便于后续复查、排错或复现：

```text
output/logs/run_YYYYMMDD_HHMMSS.log
```

## 如何导入 Web of Science 数据

1. 在 Web of Science 核心合集中使用 Topic 字段检索：

```text
TS=(false localizing sign*) OR TS=(false localising sign*) OR TS=(pseudo-localization sign)
```

2. 导出为“其他文件格式” -> “纯文本文件”或 CiteSpace 可识别的格式，字段建议包含：
   全文记录与参考文献（若做共被引分析可选）。
3. 将导出文件重命名为 `wos_citespace.txt` 并放入 `input/` 目录（或与 `config.yaml` 中 `input.wos_file` 路径一致）。
4. 重新运行 `python 0_run_all.py`。

## 输出说明

| 文件 | 说明 |
|------|------|
| `output/01_merged_records.csv` | 合并去重后的完整题录数据 |
| `output/01_citespace_input.txt` | 供 CiteSpace 读取的 WoS 格式纯文本 |
| `output/01_deduplication_report.txt` | 去重报告 |
| `output/figures/fig*.png` | 英文图表（年度趋势、期刊、国家、作者、关键词网络、词云） |
| `output/tables/*.csv` | 各类统计表格 |
| `output/report/01_False_Localizing_Sign_文献计量分析报告.docx` | 中文 Word 报告草稿（文件名前缀 01、02… 按生成先后编号，数字越大越新） |

## 注意事项

- 当前示例数据主要来自 PubMed；WoS 数据导入后将自动与之去重并合并。
- 国家/地区提取基于作者地址文本解析，可能存在少量误差，可在 `config.yaml` 的 `country_extraction.aliases` 中补充别名。
- 关键词过滤规则可在 `config.yaml` 的 `keyword.stopwords` 中调整。
- 报告中的作者、单位等信息请在 Word 中进一步补充完善。

## 依赖版本

详见 `requirements.txt`：

- pandas
- matplotlib / seaborn
- python-docx
- networkx
- wordcloud
- pycountry
- pyyaml
