#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_b1_own_data.py
用本项目自己的假定位体征文献数据（output/01_merged_records.csv）构建关键词共现网络，
生成 CiteSpace 风格的 Fig. B1（共现网络图）与 Fig. B2（聚类图），
并输出统计指标（节点/边/密度/频次/介数中心性/聚类 Q/S/LLR 标签）到 output/logs/b1_own_stats.json。
"""

import importlib.util
import itertools
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent

# ---- 复用 3_keyword_analysis.py 中的关键词清洗/映射逻辑 ----
spec = importlib.util.spec_from_file_location("kwa", ROOT / "3_keyword_analysis.py")
kwa = importlib.util.module_from_spec(spec)
spec.loader.exec_module(kwa)

config = kwa.load_config(str(ROOT / "config.yaml"))
stopwords = set(w.lower() for w in config.get("keyword", {}).get("stopwords", []))

df = pd.read_csv(ROOT / "output/01_merged_records.csv")
mapper = kwa.build_keyword_mapper(df, stopwords, config)

# 每条记录的关键词集合（去重后）
record_kws = []
for _, row in df.iterrows():
    kws = set(kwa.extract_record_keywords(row.get("keywords"), stopwords, mapper))
    record_kws.append(kws)

counts = Counter()
for kws in record_kws:
    counts.update(kws)

# ---------- 全网络（所有关键词、所有共现对）：用于报告节点/边/密度 ----------
edge_full = Counter()
for kws in record_kws:
    for a, b in itertools.combinations(sorted(kws), 2):
        edge_full[(a, b)] += 1

G_full = nx.Graph()
G_full.add_nodes_from(counts.keys())
G_full.add_edges_from(edge_full.keys())
n_nodes = G_full.number_of_nodes()
n_edges = G_full.number_of_edges()
density = nx.density(G_full)

# 介数中心性（归一化）
bc = nx.betweenness_centrality(G_full, normalized=True)

# ---------- 可视化网络：频次 >= min_count，边权 >= min_cooc ----------
min_count = config["analysis"]["min_keyword_count"]      # 3
min_cooc = config["analysis"]["min_cooccurrence"]        # 2
vis_nodes = {kw for kw, c in counts.items() if c >= min_count}
G = nx.Graph()
for kw in sorted(vis_nodes):  # 固定节点顺序，保证布局/聚类可复现
    G.add_node(kw, freq=counts[kw])
for (a, b), w in edge_full.items():
    if a in vis_nodes and b in vis_nodes and w >= min_cooc:
        G.add_edge(a, b, weight=w)
# 只保留最大连通分量用于绘图与聚类
largest = max(nx.connected_components(G), key=len)
G = G.subgraph(largest).copy()

# ---------- 聚类（greedy modularity）与 Q ----------
communities = list(nx.algorithms.community.greedy_modularity_communities(G, weight="weight"))
communities = sorted(communities, key=len, reverse=True)
comm_idx = {}
for i, com in enumerate(communities):
    for n in com:
        comm_idx[n] = i
Q = nx.algorithms.community.modularity(G, communities, weight="weight")

# ---------- 聚类标签：LLR（对数似然比） ----------
def llr(k11, k12, k21, k22):
    def xlogx(x):
        return x * math.log(x) if x > 0 else 0.0
    row1, row2 = k11 + k12, k21 + k22
    col1, col2 = k11 + k21, k12 + k22
    total = row1 + row2
    if total == 0:
        return 0.0
    e11, e12 = row1 * col1 / total, row1 * col2 / total
    e21, e22 = row2 * col1 / total, row2 * col2 / total
    return 2 * (xlogx(k11) + xlogx(k12) + xlogx(k21) + xlogx(k22)
                - xlogx(e11) - xlogx(e12) - xlogx(e21) - xlogx(e22))

# 每个簇覆盖的记录集合
cluster_records = []
for com in communities:
    recs = set()
    for i, kws in enumerate(record_kws):
        if kws & com:
            recs.add(i)
    cluster_records.append(recs)

n_records = len(record_kws)
cluster_labels = []
for ci, com in enumerate(communities):
    best, best_score = None, -1.0
    in_recs = cluster_records[ci]
    for term in com:
        term_recs = {i for i, kws in enumerate(record_kws) if term in kws}
        k11 = len(term_recs & in_recs)
        k12 = len(term_recs - in_recs)
        k21 = len(in_recs - term_recs)
        k22 = n_records - k11 - k12 - k21
        score = llr(k11, k12, k21, k22)
        if score > best_score:
            best, best_score = term, score
    cluster_labels.append(best)

# ---------- Silhouette（基于共现向量余弦距离） ----------
nodes = list(G.nodes())
adj = nx.to_numpy_array(G, nodelist=nodes, weight="weight")
norm = np.linalg.norm(adj, axis=1, keepdims=True)
norm[norm == 0] = 1.0
sim = (adj / norm) @ (adj / norm).T
dist = 1.0 - sim
np.fill_diagonal(dist, 0.0)
labels_arr = np.array([comm_idx[n] for n in nodes])

sil_per_node = np.full(len(nodes), np.nan)
for i in range(len(nodes)):
    same = labels_arr == labels_arr[i]
    same[i] = False
    if same.sum() == 0:
        continue
    a = dist[i, same].mean()
    b = min(dist[i, labels_arr == c].mean()
            for c in set(labels_arr) if c != labels_arr[i] and (labels_arr == c).sum() > 0)
    sil_per_node[i] = (b - a) / max(a, b) if max(a, b) > 0 else 0.0

cluster_sil = []
for ci in range(len(communities)):
    vals = sil_per_node[labels_arr == ci]
    vals = vals[~np.isnan(vals)]
    if len(vals):
        cluster_sil.append(float(vals.mean()))
S = float(np.mean(cluster_sil)) if cluster_sil else 0.0

# ---------- 布局（两张图共用同一坐标） ----------
pos = nx.spring_layout(G, k=1.9, iterations=400, seed=42)

freqs = np.array([G.nodes[n]["freq"] for n in G.nodes()])
bc_vis = nx.betweenness_centrality(G, normalized=True)

# ================= Fig. B1：共现网络图（CiteSpace 风格） =================
fig, ax = plt.subplots(figsize=(10, 8.4))
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

edge_ws = np.array([G.edges[e]["weight"] for e in G.edges()])
nx.draw_networkx_edges(G, pos, ax=ax, width=0.4 + 1.2 * edge_ws / edge_ws.max(),
                       edge_color="#7FB3D5", alpha=0.45)

sizes = 120 + 90 * freqs
colors = plt.cm.YlOrRd(0.35 + 0.55 * (freqs - freqs.min()) / max(1, np.ptp(freqs)))
edgecols = ["#7B2D8E" if bc_vis.get(n, 0) > 0.1 else "#C0392B" for n in G.nodes()]
nx.draw_networkx_nodes(G, pos, ax=ax, node_size=sizes, node_color=colors,
                       edgecolors=edgecols,
                       linewidths=[2.2 if bc_vis.get(n, 0) > 0.1 else 0.8 for n in G.nodes()],
                       alpha=0.9)

fmax = freqs.max()
texts = []
for n in G.nodes():
    f = G.nodes[n]["freq"]
    fs = 6.5 + 6.5 * (f / fmax)
    weight = "bold" if f >= np.percentile(freqs, 60) else "normal"
    t = ax.text(pos[n][0], pos[n][1] + 0.035 + 0.02 * (f / fmax), n,
                fontsize=fs, fontweight=weight,
                ha="center", va="center", color="#1a1a1a", zorder=5)
    texts.append(t)
kwa.adjust_labels(ax, texts)  # 防标签重叠

ax.set_title("Keyword co-occurrence network (false localizing sign research)",
             fontsize=12, pad=10)
ax.axis("off")
fig.tight_layout()
fig_b1 = ROOT / "output/figures/figB1_own_keyword_network.png"
fig.savefig(fig_b1, dpi=300, facecolor="white", bbox_inches="tight")
plt.close(fig)

# ================= Fig. B2：聚类图（彩色区块 + #N 标签） =================
def convex_hull(points):
    pts = sorted(set(map(tuple, points)))
    if len(pts) <= 2:
        return pts
    def cross(o, a, b):
        return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])
    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]

def smooth_polygon(xy, expand=0.12, subdiv=8, passes=2):
    xy = np.asarray(xy, dtype=float)
    c = xy.mean(axis=0)
    xy = c + (xy - c) * (1 + expand)
    for _ in range(subdiv):
        new = []
        n = len(xy)
        for i in range(n):
            p, q = xy[i], xy[(i + 1) % n]
            new.append(0.75 * p + 0.25 * q)
            new.append(0.25 * p + 0.75 * q)
        xy = np.array(new)
    for _ in range(passes):
        n = len(xy)
        xy = np.array([(xy[(i - 1) % n] + 2 * xy[i] + xy[(i + 1) % n]) / 4
                       for i in range(n)])
    return xy

pastel = ["#F5B7B1", "#D7BDE2", "#A9DFBF", "#F9E79F", "#AED6F1",
          "#F8C471", "#D5F5E3", "#E8DAEF", "#FCF3CF", "#D6EAF8",
          "#FADBD8", "#ABEBC6", "#FDEBD0", "#EAECEE"]
dot_colors = ["#C0392B", "#7D3C98", "#1E8449", "#B7950B", "#21618C",
              "#A04000", "#117864", "#6C3483", "#9A7D0A", "#1A5276",
              "#922B21", "#186A3B", "#7E5109", "#515A5A"]

fig, ax = plt.subplots(figsize=(10, 7.5))
fig.patch.set_facecolor("white")
ax.set_facecolor("white")

for ci, com in enumerate(communities):
    pts = np.array([pos[n] for n in com])
    if len(pts) >= 3:
        hull = convex_hull(pts)
        if len(hull) >= 3:
            poly = smooth_polygon(hull, expand=0.35)
            ax.fill(poly[:, 0], poly[:, 1], color=pastel[ci % len(pastel)],
                    alpha=0.55, zorder=1, linewidth=0)
    elif len(pts) > 0:
        c = pts.mean(axis=0)
        r = 0.18
        theta = np.linspace(0, 2 * np.pi, 60)
        ax.fill(c[0] + r * np.cos(theta), c[1] + r * np.sin(theta),
                color=pastel[ci % len(pastel)], alpha=0.55, zorder=1, linewidth=0)

label_texts = []
for ci, com in enumerate(communities):
    pts = np.array([pos[n] for n in com])
    ax.scatter(pts[:, 0], pts[:, 1], s=14, color=dot_colors[ci % len(dot_colors)],
               zorder=3, linewidths=0)
    top = pts[pts[:, 1].argmax()]
    t = ax.text(top[0], top[1] + 0.06, f"#{ci} {cluster_labels[ci]}",
                fontsize=10, fontweight="bold", color="#6B7A00",
                ha="center", va="bottom", zorder=4)
    label_texts.append(t)
kwa.adjust_labels(ax, label_texts)  # 簇标签防重叠

ax.text(0.01, 0.99,
        f"Modularity Q={Q:.3f}\nWeighted Mean Silhouette S={S:.3f}\n"
        f"N={G.number_of_nodes()}, E={G.number_of_edges()} (pruned view)",
        transform=ax.transAxes, fontsize=9, va="top", ha="left", color="#333333",
        family="monospace")
ax.set_title("Keyword cluster map (greedy modularity, LLR labels)", fontsize=12, pad=10)
ax.margins(0.12)
ax.axis("off")
fig.tight_layout()
fig_b2 = ROOT / "output/figures/figB2_own_keyword_clusters.png"
fig.savefig(fig_b2, dpi=300, facecolor="white", bbox_inches="tight")
plt.close(fig)

# ---------- 统计输出 ----------
top3_freq = counts.most_common(3)
bc_gt_01 = sorted([(k, round(v, 3)) for k, v in bc.items() if v > 0.1],
                  key=lambda x: -x[1])
cand = [(k, bc[k], counts[k]) for k in counts if counts[k] >= 5 and bc.get(k, 0) >= 0.05]
top3_bc = sorted(cand, key=lambda x: -x[1])[:3]

stats = {
    "records": n_records,
    "full_network": {
        "nodes": n_nodes,
        "links": n_edges,
        "density": round(density, 4),
    },
    "pruned_view": {
        "min_freq": min_count,
        "min_link_weight": min_cooc,
        "nodes": G.number_of_nodes(),
        "links": G.number_of_edges(),
    },
    "top3_frequency": [{"keyword": k, "freq": c} for k, c in top3_freq],
    "betweenness_gt_0.1": [{"keyword": k, "bc": v} for k, v in bc_gt_01],
    "top3_betweenness_freq_ge5_bc_ge0.05": [
        {"keyword": k, "bc": round(v, 3), "freq": f} for k, v, f in top3_bc],
    "clustering": {
        "method": "greedy modularity communities; LLR labels",
        "n_clusters": len(communities),
        "modularity_Q": round(float(Q), 4),
        "weighted_mean_silhouette_S": round(S, 4),        "clusters": [
            {"id": i, "label": cluster_labels[i], "size": len(communities[i]),
             "members_top": sorted(communities[i], key=lambda n: -counts[n])[:8]}
            for i in range(len(communities))
        ],
    },
    "figures": {"network": str(fig_b1), "clusters": str(fig_b2)},
}
out_json = ROOT / "output/logs/b1_own_stats.json"
out_json.write_text(json.dumps(stats, ensure_ascii=False, indent=2), encoding="utf-8")
print(json.dumps(stats, ensure_ascii=False, indent=2))
