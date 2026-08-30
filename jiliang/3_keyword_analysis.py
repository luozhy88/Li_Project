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


def norm_key(kw: str) -> str:
    """归一化键：小写、统一连字符/撇号、压缩空白。"""
    k = kw.lower()
    for a, b in [("–", "-"), ("—", "-"), ("‑", "-"), ("’", "'"), ("´", "'")]:
        k = k.replace(a, b)
    return re.sub(r"\s+", " ", k).strip()


def build_keyword_mapper(df: pd.DataFrame, stopwords: set, config: dict) -> dict:
    """构建 表面形式 -> 规范词 映射：先按 config 同义词表合并，
    再自动合并大小写/连字符变体（等效 VOSviewer 的 thesaurus.txt 预处理）。"""
    thesaurus = {}
    for canonical, variants in config.get("keyword", {}).get("thesaurus", {}).items():
        for v in variants:
            thesaurus[norm_key(v)] = canonical
    surfaces = Counter()
    for val in df["keywords"].fillna(""):
        for part in str(val).split(";"):
            kw = clean_keyword_text(part)
            if (kw and kw.lower() not in stopwords and len(kw) <= 100
                    and not re.match(r"^\d+$", kw)):
                surfaces[kw] += 1
    groups = {}
    for kw in surfaces:
        groups.setdefault(norm_key(kw), []).append(kw)
    mapping = {}
    for key, forms in groups.items():
        canonical = thesaurus.get(key)
        if canonical is None:
            # 频次最高者优先；并列时偏好标题式写法（非全大写/全小写）
            canonical = sorted(forms,
                               key=lambda f: (surfaces[f], not f.isupper(), not f.islower()),
                               reverse=True)[0]
        for f in forms:
            mapping[f] = canonical
    return mapping


def adjust_labels(ax, texts, iterations=1200):
    """简单斥力算法消除标签重叠（等效 VOSviewer 的 prevent label overlap）。"""
    fig = ax.figure
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    for _ in range(iterations):
        boxes = [t.get_window_extent(renderer=renderer) for t in texts]
        moved = False
        for i in range(len(texts)):
            for j in range(i + 1, len(texts)):
                b1, b2 = boxes[i], boxes[j]
                if not b1.overlaps(b2):
                    continue
                dx = min(b1.x1, b2.x1) - max(b1.x0, b2.x0)
                dy = min(b1.y1, b2.y1) - max(b1.y0, b2.y0)
                c1x, c1y = (b1.x0 + b1.x1) / 2, (b1.y0 + b1.y1) / 2
                c2x, c2y = (b2.x0 + b2.x1) / 2, (b2.y0 + b2.y1) / 2
                if dy <= dx:
                    shift = (dy / 2 + 2.5) * (1 if c1y >= c2y else -1)
                    for t, s in ((texts[i], shift), (texts[j], -shift)):
                        x, y = t.get_position()
                        y0 = inv.transform((0, 0))[1]
                        y1 = inv.transform((0, s))[1]
                        t.set_position((x, y + (y1 - y0)))
                else:
                    shift = (dx / 2 + 2.5) * (1 if c1x >= c2x else -1)
                    for t, s in ((texts[i], shift), (texts[j], -shift)):
                        x, y = t.get_position()
                        x0 = inv.transform((0, 0))[0]
                        x1 = inv.transform((s, 0))[0]
                        t.set_position((x + (x1 - x0), y))
                moved = True
        if not moved:
            break


def extract_record_keywords(keyword_string: str, stopwords: set, mapper: dict = None) -> list:
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
        result.append(mapper.get(kw, kw) if mapper else kw)
    return result


def top_keywords_table(df: pd.DataFrame, stopwords: set, config: dict, mapper: dict = None):
    all_kws = []
    for _, row in df.iterrows():
        all_kws.extend(extract_record_keywords(row.get("keywords"), stopwords, mapper))
    counts = Counter(all_kws)
    top_n = config["analysis"]["top_n_keywords"]
    top = counts.most_common(top_n)
    table = pd.DataFrame(top, columns=["Keyword", "Frequency"])
    table.to_csv(Path(config["output"]["tables"]) / "03_top_keywords.csv",
                 index=False, encoding="utf-8-sig")
    print(f"[3_keyword_analysis] Top {len(top)} keywords saved.")
    return counts, table


def cooccurrence_network(df: pd.DataFrame, counts: Counter, stopwords: set, config: dict,
                         mapper: dict = None):
    min_count = config["analysis"]["min_keyword_count"]
    min_cooc = config["analysis"]["min_cooccurrence"]

    selected = {kw for kw, c in counts.items() if c >= min_count}
    edge_counter = Counter()
    for _, row in df.iterrows():
        kws = [k for k in extract_record_keywords(row.get("keywords"), stopwords, mapper)
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
    edge_color = colors.get("grid", "#999999")
    text_color = colors.get("text", "#333333")

    node_weights = [G.nodes[n].get("weight", 1) for n in G.nodes()]

    # 模块度聚类（modularity-based clustering），不同社区赋不同颜色
    if G.number_of_edges():
        communities = list(nx.algorithms.community.greedy_modularity_communities(G))
    else:
        communities = [set(G.nodes())]
    palette = ["#4C78A8", "#F58518", "#54A24B", "#E45756", "#B279A2",
               "#EECA3B", "#72B7B2", "#9D755D", "#636363", "#17BECF"]
    comm_idx = {}
    for i, com in enumerate(sorted(communities, key=len, reverse=True)):
        for n in com:
            comm_idx[n] = i
    node_colors = [palette[comm_idx.get(n, 0) % len(palette)] for n in G.nodes()]

    fig, ax = plt.subplots(figsize=(14, 12))
    pos = nx.spring_layout(G, k=2.5, iterations=200, seed=42)

    node_sizes = [300 + 90 * w for w in node_weights]
    if G.edges(data=True):
        max_weight = max(d["weight"] for _, _, d in G.edges(data=True))
        edge_widths = [0.5 + 1.8 * d["weight"] / max_weight for _, _, d in G.edges(data=True)]
    else:
        edge_widths = []
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors,
                           alpha=0.92, edgecolors="white", linewidths=1.0, ax=ax)
    nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=0.35,
                           edge_color=edge_color, ax=ax)
    texts = nx.draw_networkx_labels(G, pos, font_size=8, font_family="sans-serif",
                                    font_color=text_color, ax=ax)
    adjust_labels(ax, list(texts.values()))  # 防标签重叠
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
    mapper = build_keyword_mapper(df, stopwords, config)
    counts, _ = top_keywords_table(df, stopwords, config, mapper)
    cooccurrence_network(df, counts, stopwords, config, mapper)
    keyword_wordcloud(counts, config)

    print("[3_keyword_analysis] Done.")


if __name__ == "__main__":
    main()
