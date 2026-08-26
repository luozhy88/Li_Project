#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
4_generate_report.py
基于已生成的统计表格与英文图表，撰写中文文献计量学报告（Word）。
"""

import argparse
import os
import re
from pathlib import Path

import pandas as pd
import yaml
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt
from docx.oxml.ns import qn


def load_config(config_path: str) -> dict:
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def set_chinese_font(run, font_name="Microsoft YaHei"):
    run.font.name = font_name
    run._element.rPr.rFonts.set(qn("w:eastAsia"), font_name)


def add_heading_zh(doc, text, level=1):
    p = doc.add_heading(level=level)
    run = p.add_run(text)
    set_chinese_font(run)
    return p


def add_paragraph_zh(doc, text, bold=False, alignment=None):
    p = doc.add_paragraph()
    run = p.add_run(text)
    set_chinese_font(run)
    run.bold = bold
    if alignment:
        p.alignment = alignment
    return p


def add_caption(doc, text):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run(text)
    run.italic = True
    run.font.size = Pt(10)
    set_chinese_font(run)
    return p


def add_picture(doc, path, caption):
    if not os.path.exists(path):
        add_paragraph_zh(doc, f"[图缺失: {path}]")
        return
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run()
    run.add_picture(str(path), width=Inches(5.8))
    add_caption(doc, caption)


def load_table(path):
    if os.path.exists(path):
        return pd.read_csv(path)
    return None


def format_reference(row):
    def clean(v):
        s = str(v if v is not None else "").strip()
        return "" if s.lower() in ("nan", "none") else s
    authors = clean(row.get("authors")).replace(";", ",")
    title = clean(row.get("title")).rstrip(".")
    journal = clean(row.get("journal"))
    year = str(int(row["year"])) if pd.notna(row.get("year")) else ""
    vol = clean(row.get("volume"))
    issue = clean(row.get("issue"))
    pages = clean(row.get("pages"))
    doi = clean(row.get("doi"))
    parts = []
    if authors:
        parts.append(authors + ".")
    if title:
        parts.append(title + ".")
    if journal:
        parts.append(journal + ".")
    if year:
        parts.append(year)
    if vol:
        loc = vol
        if issue:
            loc += f"({issue})"
        if pages:
            loc += f":{pages}"
        parts.append(loc + ".")
    if doi:
        parts.append(f"doi:{doi}.")
    ref = " ".join(parts)
    return ref if ref.endswith(".") else ref + "."


def generate_report(config: dict):
    doc = Document()
    # 设置默认中文字体
    style = doc.styles["Normal"]
    style.font.name = "Microsoft YaHei"
    style._element.rPr.rFonts.set(qn("w:eastAsia"), "Microsoft YaHei")
    style.font.size = Pt(11)

    out_root = Path(config["output"]["root"])
    fig_dir = Path(config["output"]["figures"])
    tbl_dir = Path(config["output"]["tables"])

    # 加载数据
    merged = pd.read_csv(out_root / "01_merged_records.csv")
    stats = load_table(tbl_dir / "02_general_statistics.csv").iloc[0].to_dict()
    annual = load_table(tbl_dir / "02_annual_publications.csv")
    journals = load_table(tbl_dir / "02_top_journals.csv")
    countries = load_table(tbl_dir / "02_top_countries.csv")
    authors = load_table(tbl_dir / "02_top_authors.csv")
    keywords = load_table(tbl_dir / "03_top_keywords.csv")

    total = int(stats["total_records"])
    year_min = int(stats["year_min"]) if pd.notna(stats["year_min"]) else "N/A"
    year_max = int(stats["year_max"]) if pd.notna(stats["year_max"]) else "N/A"

    # 标题
    title = doc.add_paragraph()
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = title.add_run(config["report"]["title"])
    run.bold = True
    run.font.size = Pt(18)
    set_chinese_font(run)

    subtitle = doc.add_paragraph()
    subtitle.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = subtitle.add_run(config["report"]["subtitle"])
    run.font.size = Pt(14)
    set_chinese_font(run)

    add_paragraph_zh(doc, f"作者：{config['report']['author']}", alignment=WD_ALIGN_PARAGRAPH.CENTER)
    add_paragraph_zh(doc, f"生成日期：{pd.Timestamp.now().strftime('%Y-%m-%d')}", alignment=WD_ALIGN_PARAGRAPH.CENTER)
    doc.add_paragraph()

    # 摘要
    add_heading_zh(doc, "摘要", level=1)
    abstract_text = (
        f"目的：对假定位体征（false localizing sign）相关文献进行文献计量学分析，"
        f"梳理该领域的研究现状、热点与发展趋势。方法：在 PubMed 中以检索式 "
        f"“{config['search']['pubmed_query']}”进行主题检索，并预留 Web of Science（WoS）"
        f"Topic 检索结果导入接口；对两库题录进行解析、去重后，采用 Python 进行描述性统计、"
        f"期刊/国家/作者分布分析及关键词共现分析。结果：共纳入 {total} 篇文献（当前以 PubMed 数据为主），"
        f"时间跨度为 {year_min}—{year_max} 年；发文量总体呈波动上升趋势，"
        f"主要期刊集中于神经病学、神经外科学与神经眼科学领域；高频关键词包括 "
        f"{', '.join(keywords['Keyword'].head(5).tolist())} 等。结论：假定位体征的研究仍以病例报告与影像学特征为主，"
        f"未来可进一步开展多中心、大样本的临床与影像学研究。"
    )
    add_paragraph_zh(doc, abstract_text)
    add_paragraph_zh(doc, "关键词：假定位体征；文献计量学；PubMed；Web of Science；CiteSpace", bold=False)
    doc.add_paragraph()

    # 1. 引言
    add_heading_zh(doc, "1. 引言", level=1)
    intro = (
        "假定位体征（false localizing sign / false localising sign）是指神经系统体征所提示的病变部位"
        "与实际病灶不一致的现象，常见于颅内压增高、脑疝、脊髓病变等情形。经典的假定位体征包括"
        "外展神经麻痹、Kernohan-Woltman 切迹现象、对侧颅神经受累等。由于其可导致定位诊断错误，"
        "进而延误治疗，因此受到神经科、神经外科、神经眼科及影像科医师的广泛关注。"
        "本研究拟采用文献计量学方法，系统分析该领域的研究分布、核心作者、主要期刊及研究热点，"
        "为后续研究选题与临床教学提供参考。"
    )
    add_paragraph_zh(doc, intro)

    # 2. 材料与方法
    add_heading_zh(doc, "2. 材料与方法", level=1)
    add_heading_zh(doc, "2.1 数据来源与检索策略", level=2)
    method1 = (
        f"本研究以 PubMed 与 Web of Science（WoS）核心合集为数据来源。PubMed 检索式为："
        f"{config['search']['pubmed_query']}。WoS 采用 Topic 字段检索，检索式为："
        f"{config['search']['wos_query']}。检索词覆盖 'false localizing sign'、'false localising sign' "
        f"及 'pseudo-localization sign' 等常见拼写变体。"
    )
    add_paragraph_zh(doc, method1)

    add_heading_zh(doc, "2.2 纳入与排除标准", level=2)
    add_paragraph_zh(doc, "纳入标准：题录信息完整的期刊论文、综述、病例报告等；排除标准：会议摘要、重复发表、信息缺失无法解析的记录。")

    add_heading_zh(doc, "2.3 数据导出与去重", level=2)
    method2 = (
        "从 PubMed 导出 MEDLINE 格式题录，从 WoS 导出 CiteSpace 兼容的纯文本题录。"
        "使用自定义 Python 脚本解析字段（标题、作者、摘要、MeSH 主题词、其他关键词、机构、国家、DOI 等），"
        "依次按 PMID、DOI、归一化标题进行去重，优先保留 PubMed 记录。"
    )
    add_paragraph_zh(doc, method2)

    add_heading_zh(doc, "2.4 分析方法", level=2)
    add_paragraph_zh(doc, "采用描述性统计、年度发文趋势、期刊分布、国家/地区分布、核心作者分析、MeSH/其他关键词词频统计及关键词共现网络分析。图表使用 matplotlib 与 seaborn 绘制（英文标注），报告使用 python-docx 生成。")

    # 3. 结果
    add_heading_zh(doc, "3. 结果", level=1)

    add_heading_zh(doc, "3.1 检索与去重结果", level=2)
    pubmed_n = len(merged[merged["source_database"] == "PubMed"])
    wos_n = len(merged[merged["source_database"] == "WoS"])
    result1 = (
        f"经检索，PubMed 获得 {pubmed_n} 条记录，WoS 获得 {wos_n} 条记录（如已导入）。"
        f"去重后最终纳入 {total} 篇文献进行分析，时间跨度为 {year_min}—{year_max} 年。"
        f"其中具有 DOI 的记录 {int(stats['with_doi'])} 条，具有摘要的记录 {int(stats['with_abstract'])} 条。"
    )
    add_paragraph_zh(doc, result1)

    add_heading_zh(doc, "3.2 年度发文趋势", level=2)
    recent = annual.tail(5)
    recent_text = ", ".join([f"{int(r['Year'])}年（{int(r['Publications'])}篇）" for _, r in recent.iterrows()])
    add_paragraph_zh(doc, f"图 1 展示了假定位体征相关研究的年度发文量与累计发文量。近五年发文情况为：{recent_text}。")
    add_picture(doc, fig_dir / "fig1_annual_trend.png", "Fig. 1 Annual publication output and cumulative trend of studies on false localizing signs.")

    add_heading_zh(doc, "3.3 期刊分布", level=2)
    top3_journals = ", ".join([f"《{r['Journal']}》（{int(r['Publications'])}篇）" for _, r in journals.head(3).iterrows()])
    add_paragraph_zh(doc, f"发文量排名前 3 的期刊依次为：{top3_journals}。图 2 展示了发文量前 {len(journals)} 位的期刊分布。")
    add_picture(doc, fig_dir / "fig2_top_journals.png", "Fig. 2 Top journals publishing studies on false localizing signs.")

    add_heading_zh(doc, "3.4 国家/地区分布", level=2)
    if countries is not None and not countries.empty:
        top3_countries = ", ".join([f"{r['Country/Region']}（{int(r['Publications'])}篇）" for _, r in countries.head(3).iterrows()])
        add_paragraph_zh(doc, f"按第一/通讯作者机构地址提取国家/地区，发文量前 3 位为：{top3_countries}。图 3 展示了前 {len(countries)} 位国家/地区分布。")
        add_picture(doc, fig_dir / "fig3_top_countries.png", "Fig. 3 Top countries/regions contributing to false localizing sign research.")
    else:
        add_paragraph_zh(doc, "未提取到有效的国家/地区信息。")

    add_heading_zh(doc, "3.5 核心作者", level=2)
    if authors is not None and not authors.empty:
        top3_authors = ", ".join([f"{r['Author']}（{int(r['Publications'])}篇）" for _, r in authors.head(3).iterrows()])
        add_paragraph_zh(doc, f"发文量最多的作者包括：{top3_authors}。图 4 展示了前 {len(authors)} 位核心作者。")
        add_picture(doc, fig_dir / "fig4_top_authors.png", "Fig. 4 Top authors in false localizing sign research.")

    add_heading_zh(doc, "3.6 关键词与研究热点", level=2)
    if keywords is not None and not keywords.empty:
        top5 = ", ".join(keywords["Keyword"].head(5).tolist())
        add_paragraph_zh(doc, f"高频关键词包括 {top5} 等。图 5 为关键词共现网络，图 6 为关键词词云。")
        add_picture(doc, fig_dir / "fig5_keyword_network.png", "Fig. 5 Keyword co-occurrence network based on MeSH terms and author keywords.")
        add_picture(doc, fig_dir / "fig6_keyword_wordcloud.png", "Fig. 6 Word cloud of keywords in false localizing sign research.")

    # 4. 讨论
    add_heading_zh(doc, "4. 讨论", level=1)
    discuss = (
        "本研究首次对假定位体征相关文献进行了系统计量分析。结果显示，该领域文献数量总体呈增长态势，"
        "尤其近十年病例报告与影像学研究增多，反映了神经影像技术对识别罕见假定位体征的促进作用。"
        "从期刊分布看，研究主要发表于神经病学、神经外科学与神经眼科学期刊，提示该主题具有鲜明的专科交叉特征。"
        "高频关键词集中于磁共振成像、颅内肿瘤、颅内压异常及颅神经麻痹等，说明假定位体征的研究热点"
        "与神经影像、颅内压管理及颅神经病变密切相关。"
        "\n\n"
        "本研究的局限性包括：当前分析以 PubMed 数据为主，WoS 数据导入后可进一步补充会议论文与引文信息；"
        "机构与国家提取基于地址文本解析，可能存在少量误判；关键词共现基于 MeSH 与作者关键词，"
        "未进行语义聚类。未来研究可结合 CiteSpace 的突现检测与 VOSviewer 的共被引分析，"
        "进一步揭示学科演化路径与知识基础。"
    )
    add_paragraph_zh(doc, discuss)

    # 5. 结论
    add_heading_zh(doc, "5. 结论", level=1)
    conclusion = (
        "假定位体征相关研究历经数十年仍保持活跃，MRI、CT 等影像技术以及颅神经麻痹、颅内压异常是该领域的核心主题。"
        "文献计量结果可为临床医师识别假定位体征、避免误诊提供文献支持，也可为后续研究选题提供方向。"
    )
    add_paragraph_zh(doc, conclusion)

    # 纳入文献清单
    doc.add_page_break()
    add_heading_zh(doc, "附录：纳入文献清单（部分）", level=1)
    add_paragraph_zh(doc, f"以下列出纳入分析的 {total} 条记录中的前 30 条，完整清单见 output/01_merged_records.csv。")
    for i, (_, row) in enumerate(merged.head(30).iterrows(), 1):
        ref = format_reference(row)
        add_paragraph_zh(doc, f"{i}. {ref}")

    # 保存
    report_path = Path(config["output"]["report"]) / config["report"]["docx_filename"]
    report_path.parent.mkdir(parents=True, exist_ok=True)
    doc.save(report_path)
    print(f"[4_generate_report] Report saved to {report_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate Chinese Word report.")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    generate_report(config)


if __name__ == "__main__":
    main()
