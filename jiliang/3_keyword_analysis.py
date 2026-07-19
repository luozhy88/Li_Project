#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
3_keyword_analysis.py
关键词分析：高频 MeSH/OT 关键词、关键词共现网络、关键词词云。
"""

import argparse
import itertools
import re
from collections import Counter
from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import seaborn as sns
import yaml
from wordcloud import WordCloud


def load_config(config_path: str) -> dict:
    with open(config_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def ensure_dirs(config: dict):
    for d in [config["output"]["root"], config["output"]["figures"],
              config["output"]["tables"]]:
        Path(d).mkdir(parents=True, exist_ok=True)


def setup_matplotlib_style(config: dict):
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
        "figure.dpi": config["analysis"]["figure_dpi"],
        "savefig.dpi": config["analysis"]["figure_dpi"],
        "savefig.facecolor": "white",
        "savefig.bbox": "tight",
    })


def clean_keyword_text(kw: str) -> str:
    if not kw:
        return ""
    kw = kw.strip()
    # 去掉末尾的 * 标记
    kw = kw.rstrip("*").strip()
    return kw


def extract_record_keywords(keyword_string: str, stopwords: set) -> list:
    if pd.isna(keyword_string) or not str(keyword_string).strip():
        return []
    result = []
    for part in str(keyword_string).split(";"):
        kw = clean_keyword_text(part)
        if not kw:
            continue
        if kw.lower() in stopwords:
            continue
        # 排除过长或纯数字
        if len(kw) > 100 or re.match(r"^\d+$", kw):
            continue
        result.append(kw)
    return result


def top_keywords_table(df: pd.DataFrame, stopwords: set, config: dict):
    all_kws = []
    for _, row in df.iterrows():
        all_kws.extend(extract_record_keywords(row.get("keywords"), stopwords))
    counts = Counter(all_kws)
    top_n = config["analysis"]["top_n_keywords"]
    top = counts.most_common(top_n)
    table = pd.DataFrame(top, columns=["Keyword", "Frequency"])
    table.to_csv(Path(config["output"]["tables"]) / "03_top_keywords.csv",
                 index=False, encoding="utf-8-sig")
    print(f"[3_keyword_analysis] Top {len(top)} keywords saved.")
    return counts, table


def cooccurrence_network(df: pd.DataFrame, counts: Counter, stopwords: set, config: dict):
    min_count = config["analysis"]["min_keyword_count"]
    min_cooc = config["analysis"]["min_cooccurrence"]

    selected = {kw for kw, c in counts.items() if c >= min_count}
    edge_counter = Counter()
    for _, row in df.iterrows():
        kws = [k for k in extract_record_keywords(row.get("keywords"), stopwords)
               if k in selected]
        for a, b in itertools.combinations(sorted(set(kws)), 2):
            edge_counter[(a, b)] += 1

    edges = [(a, b, w) for (a, b), w in edge_counter.items() if w >= min_cooc]
    edges_df = pd.DataFrame(edges, columns=["Source", "Target", "Weight"])
    edges_df.to_csv(Path(config["output"]["tables"]) / "03_keyword_edges.csv",
                    index=False, encoding="utf-8-sig")

    node_freq = {kw: counts[kw] for kw in selected if counts[kw] >= min_count}
    nodes_df = pd.DataFrame(list(node_freq.items()), columns=["Keyword", "Frequency"])
    nodes_df.to_csv(Path(config["output"]["tables"]) / "03_keyword_nodes.csv",
                    index=False, encoding="utf-8-sig")

    # 绘图
    G = nx.Graph()
    for kw, freq in node_freq.items():
        G.add_node(kw, weight=freq)
    for a, b, w in edges:
        G.add_edge(a, b, weight=w)

    colors = config.get("colors", {})
    primary = colors.get("primary", "#2874A6")
    edge_color = colors.get("grid", "#999999")
    text_color = colors.get("text", "#333333")

    # 节点按频次使用从浅到深的渐变（低频次仍保持可见）
    node_weights = [G.nodes[n].get("weight", 1) for n in G.nodes()]
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "nature_nodes", ["#D6EAF8", primary])

    fig, ax = plt.subplots(figsize=(14, 12))
    pos = nx.spring_layout(G, k=2.5, iterations=200, seed=42)

    node_sizes = [300 + 90 * w for w in node_weights]
    if G.edges(data=True):
        max_weight = max(d["weight"] for _, _, d in G.edges(data=True))
        edge_widths = [0.5 + 1.8 * d["weight"] / max_weight for _, _, d in G.edges(data=True)]
    else:
        edge_widths = []
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_weights,
                           cmap=cmap, alpha=0.92, edgecolors="white", linewidths=1.0, ax=ax)
    nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.35,
                           edge_color=edge_color, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=8, font_family="sans-serif",
                            font_color=text_color, ax=ax)
    ax.set_title("Keyword co-occurrence network (MeSH / other terms)")
    ax.axis("off")
    fig.tight_layout()
    out_path = Path(config["output"]["figures"]) / "fig5_keyword_network.png"
    fig.savefig(out_path, dpi=config["analysis"]["figure_dpi"],
                format=config["analysis"]["figure_format"])
    plt.close(fig)
    print(f"[3_keyword_analysis] Saved {out_path}")


def keyword_wordcloud(counts: Counter, config: dict):
    if not counts:
        return
    colors = config.get("colors", {})
    wc_palette = colors.get("wordcloud", ["#0072B2", "#D55E00", "#009E73",
                                          "#CC79A7", "#56B4E9", "#E69F00"])
    cmap = mcolors.ListedColormap(wc_palette)
    wc = WordCloud(width=1200, height=600, background_color="white",
                   colormap=cmap, max_words=100,
                   prefer_horizontal=0.9).generate_from_frequencies(counts)
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.imshow(wc, interpolation="bilinear")
    ax.axis("off")
    ax.set_title("Word cloud of keywords")
    fig.tight_layout()
    out_path = Path(config["output"]["figures"]) / "fig6_keyword_wordcloud.png"
    fig.savefig(out_path, dpi=config["analysis"]["figure_dpi"],
                format=config["analysis"]["figure_format"])
    plt.close(fig)
    print(f"[3_keyword_analysis] Saved {out_path}")


def main():
    parser = argparse.ArgumentParser(description="Keyword analysis for bibliometric data.")
    parser.add_argument("--config", default="config.yaml")
    args = parser.parse_args()

    config = load_config(args.config)
    ensure_dirs(config)
    sns.set_style(config["analysis"]["seaborn_style"])
    setup_matplotlib_style(config)

    merged_path = Path(config["output"]["root"]) / "01_merged_records.csv"
    df = pd.read_csv(merged_path)
    print(f"[3_keyword_analysis] Loaded {len(df)} records.")

    stopwords = set([w.lower() for w in config.get("keyword", {}).get("stopwords", [])])
    counts, _ = top_keywords_table(df, stopwords, config)
    cooccurrence_network(df, counts, stopwords, config)
    keyword_wordcloud(counts, config)

    print("[3_keyword_analysis] Done.")


if __name__ == "__main__":
    main()
