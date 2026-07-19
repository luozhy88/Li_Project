#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2_descriptive_analysis.py
对合并后的文献题录进行描述性文献计量分析：
年度趋势、期刊分布、国家/地区分布、作者分布、文献类型分布。
输出表格（CSV）与英文图表（PNG）。
"""

import argparse
import os
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import pycountry
import seaborn as sns
import yaml


def load_config(config_path: str) -> dict:
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def ensure_dirs(config: dict):
    for d in [config["output"]["root"], config["output"]["figures"],
              config["output"]["tables"]]:
        Path(d).mkdir(parents=True, exist_ok=True)


def setup_matplotlib_style(config: dict):
    """应用 Nature/Science 风格的简洁图表样式。"""
    colors = config.get("colors", {})
    text_color = colors.get("text", "#333333")
    grid_color = colors.get("grid", "#E5E5E5")
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "axes.edgecolor": text_color,
        "axes.labelcolor": text_color,
        "axes.linewidth": 0.8,
        "xtick.color": text_color,
        "ytick.color": text_color,
        "text.color": text_color,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.grid": True,
        "axes.grid.axis": "y",
        "grid.color": grid_color,
        "grid.linewidth": 0.8,
        "grid.linestyle": "-",
        "figure.dpi": config["analysis"]["figure_dpi"],
        "savefig.dpi": config["analysis"]["figure_dpi"],
        "savefig.facecolor": "white",
        "savefig.bbox": "tight",
    })


def build_country_patterns(config: dict):
    """构建国家/地区名称匹配模式（含别名），长名称优先，输出统一规范名。"""
    aliases = config.get("country_extraction", {}).get("aliases", {})
    mapping = {}  # 任意名称/别名 -> 规范名
    for country in pycountry.countries:
        canonical = aliases.get(country.name, country.name)
        for n in {country.name,
                  getattr(country, "official_name", ""),
                  getattr(country, "common_name", "")}:
            if not n:
                continue
            mapping[n] = aliases.get(n, canonical)
    # 显式别名 -> 规范名
    for alias, canonical in aliases.items():
        mapping[alias] = canonical

    # 生成 (canonical, pattern) 列表，按 pattern 长度降序，避免短名提前匹配
    items = []
    seen_canonical = set()
    for name, canonical in mapping.items():
        if canonical in seen_canonical:
            continue
        items.append((canonical, re.escape(name)))
        seen_canonical.add(canonical)
    # 补充尚未加入的别名 pattern（同名 canonical 取最长 pattern）
    extra = []
    for name, canonical in mapping.items():
        extra.append((canonical, re.escape(name)))
    # 合并去重，按 pattern 长度降序
    unique = {}
    for canonical, pat in items + extra:
        if pat not in unique:
            unique[pat] = canonical
    sorted_items = sorted(
        [(canonical, pat) for pat, canonical in unique.items()],
        key=lambda x: len(x[1]),
        reverse=True,
    )
    return sorted_items


def extract_countries(affiliation: str, country_patterns) -> list:
    """从地址字符串中提取国家/地区，返回规范名列表。"""
    if not affiliation:
        return []
    text = affiliation.lower()
    found = []
    for canonical, pattern in country_patterns:
        if re.search(r"\b" + pattern.lower() + r"\b", text):
            if canonical not in found:
                found.append(canonical)
            # 替换已匹配部分，避免子串重复匹配
            text = re.sub(r"\b" + pattern.lower() + r"\b", "", text)
    return found


def plot_annual_trend(df: pd.DataFrame, config: dict):
    colors = config.get("colors", {})
    primary = colors.get("primary", "#2874A6")
    secondary = colors.get("secondary", "#C0392B")

    counts = df["year"].value_counts().sort_index()
    counts = counts[counts.index.notna()]
    counts.index = counts.index.astype(int)
    cumulative = counts.cumsum()

    fig, ax1 = plt.subplots(figsize=tuple(config["analysis"]["figsize"]))
    ax1.bar(counts.index, counts.values, color=primary, alpha=0.85,
            edgecolor="white", linewidth=0.5, label="Publications")
    ax1.set_xlabel("Year")
    ax1.set_ylabel("Number of publications")
    ax1.set_title("Annual publication output on false localizing signs")
    ax1.tick_params(axis="y")

    ax2 = ax1.twinx()
    ax2.plot(counts.index, cumulative.values, color=secondary, marker="o",
             markersize=5, linewidth=2, label="Cumulative")
    ax2.set_ylabel("Cumulative publications")
    ax2.tick_params(axis="y")

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc="upper left",
               frameon=False)
    sns.despine(ax=ax1, top=True, right=False)
    sns.despine(ax=ax2, top=True, right=False, left=True)

    fig.tight_layout()
    out_path = Path(config["output"]["figures"]) / "fig1_annual_trend.png"
    fig.savefig(out_path, dpi=config["analysis"]["figure_dpi"],
                format=config["analysis"]["figure_format"])
    plt.close(fig)
    print(f"[2_descriptive_analysis] Saved {out_path}")

    table = pd.DataFrame({"Year": counts.index, "Publications": counts.values,
                          "Cumulative": cumulative.values})
    table.to_csv(Path(config["output"]["tables"]) / "02_annual_publications.csv",
                 index=False, encoding="utf-8-sig")
    return table


def plot_top_journals(df: pd.DataFrame, config: dict):
    colors = config.get("colors", {})
    primary = colors.get("primary", "#2874A6")
    top_n = config["analysis"]["top_n_journals"]
    counts = df["journal"].value_counts().head(top_n)

    fig, ax = plt.subplots(figsize=(config["analysis"]["figsize"][0],
                                     max(6, top_n * 0.4)))
    sns.barplot(x=counts.values, y=counts.index, color=primary, ax=ax)
    ax.set_xlabel("Number of publications")
    ax.set_ylabel("Journal")
    ax.set_title(f"Top {len(counts)} journals publishing on false localizing signs")
    sns.despine(ax=ax, top=True, right=True)
    fig.tight_layout()
    out_path = Path(config["output"]["figures"]) / "fig2_top_journals.png"
    fig.savefig(out_path, dpi=config["analysis"]["figure_dpi"],
                format=config["analysis"]["figure_format"])
    plt.close(fig)
    print(f"[2_descriptive_analysis] Saved {out_path}")

    table = counts.reset_index()
    table.columns = ["Journal", "Publications"]
    table.to_csv(Path(config["output"]["tables"]) / "02_top_journals.csv",
                 index=False, encoding="utf-8-sig")
    return table


def plot_top_countries(df: pd.DataFrame, country_patterns, config: dict):
    colors = config.get("colors", {})
    primary = colors.get("primary", "#2874A6")
    top_n = config["analysis"]["top_n_countries"]
    country_records = []
    for _, row in df.iterrows():
        aff = str(row.get("affiliations", ""))
        countries = extract_countries(aff, country_patterns)
        for c in countries:
            country_records.append({"pmid": row.get("pmid"), "country": c})

    if country_records:
        country_df = pd.DataFrame(country_records)
        counts = country_df["country"].value_counts().head(top_n)
    else:
        counts = pd.Series(dtype=int)

    fig, ax = plt.subplots(figsize=tuple(config["analysis"]["figsize"]))
    if not counts.empty:
        sns.barplot(x=counts.values, y=counts.index, color=primary, ax=ax)
    ax.set_xlabel("Number of publications")
    ax.set_ylabel("Country / Region")
    ax.set_title(f"Top {len(counts)} countries/regions publishing on false localizing signs")
    sns.despine(ax=ax, top=True, right=True)
    fig.tight_layout()
    out_path = Path(config["output"]["figures"]) / "fig3_top_countries.png"
    fig.savefig(out_path, dpi=config["analysis"]["figure_dpi"],
                format=config["analysis"]["figure_format"])
    plt.close(fig)
    print(f"[2_descriptive_analysis] Saved {out_path}")

    table = counts.reset_index()
    table.columns = ["Country/Region", "Publications"]
    table.to_csv(Path(config["output"]["tables"]) / "02_top_countries.csv",
                 index=False, encoding="utf-8-sig")
    return table, country_df if country_records else pd.DataFrame()


def plot_top_authors(df: pd.DataFrame, config: dict):
    colors = config.get("colors", {})
    primary = colors.get("primary", "#2874A6")
    top_n = config["analysis"]["top_n_authors"]
    authors = []
    for val in df["authors"].dropna():
        for a in str(val).split(";"):
            a = a.strip()
            if a:
                authors.append(a)
    counts = pd.Series(authors).value_counts().head(top_n)

    fig, ax = plt.subplots(figsize=(config["analysis"]["figsize"][0],
                                     max(6, top_n * 0.4)))
    sns.barplot(x=counts.values, y=counts.index, color=primary, ax=ax)
    ax.set_xlabel("Number of publications")
    ax.set_ylabel("Author")
    ax.set_title(f"Top {len(counts)} most productive authors")
    sns.despine(ax=ax, top=True, right=True)
    fig.tight_layout()
    out_path = Path(config["output"]["figures"]) / "fig4_top_authors.png"
    fig.savefig(out_path, dpi=config["analysis"]["figure_dpi"],
                format=config["analysis"]["figure_format"])
    plt.close(fig)
    print(f"[2_descriptive_analysis] Saved {out_path}")

    table = counts.reset_index()
    table.columns = ["Author", "Publications"]
    table.to_csv(Path(config["output"]["tables"]) / "02_top_authors.csv",
                 index=False, encoding="utf-8-sig")
    return table


def publication_type_table(df: pd.DataFrame, config: dict):
    types = []
    for val in df["publication_types"].dropna():
        for t in str(val).split(";"):
            t = t.strip()
            if t:
                types.append(t)
    counts = pd.Series(types).value_counts()
    table = counts.reset_index()
    table.columns = ["Publication type", "Count"]
    table.to_csv(Path(config["output"]["tables"]) / "02_publication_types.csv",
                 index=False, encoding="utf-8-sig")
    return table


def main():
    parser = argparse.ArgumentParser(description="Descriptive bibliometric analysis.")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    ensure_dirs(config)

    sns.set_style(config["analysis"]["seaborn_style"])
    setup_matplotlib_style(config)

    merged_path = Path(config["output"]["root"]) / "01_merged_records.csv"
    if not merged_path.exists():
        print(f"[2_descriptive_analysis] ERROR: {merged_path} not found. Run 1_parse_and_merge.py first.")
        sys.exit(1)

    df = pd.read_csv(merged_path)
    print(f"[2_descriptive_analysis] Loaded {len(df)} records.")

    # 全局统计
    stats = {
        "total_records": len(df),
        "year_min": int(df["year"].min()) if df["year"].notna().any() else None,
        "year_max": int(df["year"].max()) if df["year"].notna().any() else None,
        "with_doi": df["doi"].notna().sum(),
        "with_abstract": df["abstract"].notna().sum(),
    }
    stats_df = pd.DataFrame([stats])
    stats_df.to_csv(Path(config["output"]["tables"]) / "02_general_statistics.csv",
                    index=False, encoding="utf-8-sig")
    print("[2_descriptive_analysis] General statistics:")
    for k, v in stats.items():
        print(f"  {k}: {v}")

    country_patterns = build_country_patterns(config)

    plot_annual_trend(df, config)
    plot_top_journals(df, config)
    _, _ = plot_top_countries(df, country_patterns, config)
    plot_top_authors(df, config)
    publication_type_table(df, config)

    print("[2_descriptive_analysis] Done.")


if __name__ == "__main__":
    main()
