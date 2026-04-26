#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
网络拓扑分析、聚类分析、富集分析、显著性检验
基于清理后的数据库数据
"""

import sqlite3
import json
import csv
import os
import requests
import time
from collections import defaultdict
import networkx as nx
from networkx.algorithms import community
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np

OUTPUT_DIR = 'output/analysis_results'
os.makedirs(OUTPUT_DIR, exist_ok=True)


def build_network(conn):
    """构建穴位-基因二分网络"""
    cursor = conn.cursor()
    cursor.execute('''
        SELECT a.name_cn, g.gene_symbol, COUNT(*) as cnt
        FROM acupoints a
        JOIN article_acupoint aa ON a.id = aa.acupoint_id
        JOIN article_gene ag ON aa.article_id = ag.article_id
        JOIN genes g ON ag.gene_id = g.id
        WHERE g.quality_flag IN ('valid', 'mapped_to_official')
        GROUP BY a.id, g.id
        HAVING cnt >= 2
        ORDER BY cnt DESC
    ''')
    
    G = nx.Graph()
    acupoint_nodes = set()
    gene_nodes = set()
    edges = []
    
    for row in cursor.fetchall():
        apt, gene, count = row
        G.add_edge(apt, gene, weight=count)
        acupoint_nodes.add(apt)
        gene_nodes.add(gene)
        edges.append((apt, gene, count))
    
    return G, acupoint_nodes, gene_nodes, edges


def network_topology_analysis(G, acupoint_nodes, gene_nodes):
    """网络拓扑分析：度中心性、介数中心性、接近中心性"""
    print("=" * 60)
    print("1. 网络拓扑分析")
    print("=" * 60)
    
    results = {
        'network_stats': {},
        'degree_centrality': {},
        'betweenness_centrality': {},
        'closeness_centrality': {},
        'eigenvector_centrality': {},
    }
    
    # 基本网络统计
    results['network_stats'] = {
        'total_nodes': G.number_of_nodes(),
        'total_edges': G.number_of_edges(),
        'density': round(nx.density(G), 4),
        'average_clustering': round(nx.average_clustering(G), 4),
        'connected_components': nx.number_connected_components(G),
    }
    
    print(f"网络节点数: {results['network_stats']['total_nodes']}")
    print(f"网络边数: {results['network_stats']['total_edges']}")
    print(f"网络密度: {results['network_stats']['density']}")
    print(f"平均聚类系数: {results['network_stats']['average_clustering']}")
    print(f"连通分量数: {results['network_stats']['connected_components']}")
    
    # 度中心性（加权）
    degree_cent = nx.degree_centrality(G)
    weighted_degree = dict(G.degree(weight='weight'))
    
    # 介数中心性
    betweenness_cent = nx.betweenness_centrality(G, weight='weight')
    
    # 接近中心性
    closeness_cent = nx.closeness_centrality(G)
    
    # 特征向量中心性
    try:
        eigenvector_cent = nx.eigenvector_centrality(G, weight='weight', max_iter=1000)
    except:
        eigenvector_cent = {}
    
    # 分别统计穴位和基因的中心性
    apt_degree = {n: degree_cent[n] for n in acupoint_nodes}
    gene_degree = {n: degree_cent[n] for n in gene_nodes}
    
    apt_betweenness = {n: betweenness_cent[n] for n in acupoint_nodes}
    gene_betweenness = {n: betweenness_cent[n] for n in gene_nodes}
    
    results['degree_centrality']['top_acupoints'] = sorted(apt_degree.items(), key=lambda x: -x[1])[:10]
    results['degree_centrality']['top_genes'] = sorted(gene_degree.items(), key=lambda x: -x[1])[:15]
    results['betweenness_centrality']['top_acupoints'] = sorted(apt_betweenness.items(), key=lambda x: -x[1])[:10]
    results['betweenness_centrality']['top_genes'] = sorted(gene_betweenness.items(), key=lambda x: -x[1])[:15]
    results['weighted_degree'] = sorted(weighted_degree.items(), key=lambda x: -x[1])[:20]
    
    print("\n--- 度中心性 Top 穴位 ---")
    for node, cent in results['degree_centrality']['top_acupoints']:
        print(f"  {node}: {cent:.4f}")
    
    print("\n--- 度中心性 Top 基因 ---")
    for node, cent in results['degree_centrality']['top_genes']:
        print(f"  {node}: {cent:.4f}")
    
    print("\n--- 介数中心性 Top 基因 ---")
    for node, cent in results['betweenness_centrality']['top_genes']:
        print(f"  {node}: {cent:.4f}")
    
    print("\n--- 加权度 Top 节点 ---")
    for node, deg in results['weighted_degree'][:15]:
        print(f"  {node}: {deg}")
    
    return results


def cluster_analysis(G, acupoint_nodes, gene_nodes):
    """聚类分析：Louvain社区发现"""
    print("\n" + "=" * 60)
    print("2. 聚类分析 (Louvain Community Detection)")
    print("=" * 60)
    
    try:
        import community as community_louvain
        partition = community_louvain.best_partition(G, weight='weight')
    except ImportError:
        print("python-louvain not available, using greedy modularity communities")
        communities = community.greedy_modularity_communities(G, weight='weight')
        partition = {}
        for i, comm in enumerate(communities):
            for node in comm:
                partition[node] = i
    
    # 统计每个社区
    communities = defaultdict(list)
    for node, comm_id in partition.items():
        communities[comm_id].append(node)
    
    results = {
        'num_communities': len(communities),
        'communities': []
    }
    
    print(f"检测到 {len(communities)} 个社区")
    
    for comm_id, nodes in sorted(communities.items(), key=lambda x: -len(x[1])):
        apt_in_comm = [n for n in nodes if n in acupoint_nodes]
        genes_in_comm = [n for n in nodes if n in gene_nodes]
        
        comm_info = {
            'community_id': comm_id,
            'size': len(nodes),
            'acupoints': apt_in_comm,
            'genes': genes_in_comm,
        }
        results['communities'].append(comm_info)
        
        if len(apt_in_comm) > 0:
            print(f"\n社区 {comm_id}: {len(nodes)} 个节点")
            print(f"  穴位: {', '.join(apt_in_comm)}")
            print(f"  基因: {', '.join(genes_in_comm[:10])}{'...' if len(genes_in_comm) > 10 else ''}")
    
    return results, partition


def enrichment_analysis(gene_list, description=""):
    """使用Enrichr API进行GO/KEGG富集分析"""
    print("\n" + "=" * 60)
    print("3. GO/KEGG 富集分析 (Enrichr API)")
    print("=" * 60)
    
    if not gene_list:
        print("基因列表为空，跳过富集分析")
        return None
    
    # Enrichr API endpoint
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr'
    
    results = {}
    
    try:
        # 添加基因列表到Enrichr
        add_list_url = f"{ENRICHR_URL}/addList"
        genes_str = '\n'.join(gene_list)
        response = requests.post(add_list_url, files={
            'list': (None, genes_str),
            'description': (None, description)
        }, timeout=30)
        
        if response.status_code != 200:
            print(f"Enrichr addList failed: {response.status_code}")
            return None
        
        data = response.json()
        user_list_id = data.get('userListId')
        
        if not user_list_id:
            print("No userListId returned from Enrichr")
            return None
        
        print(f"Enrichr list added, userListId: {user_list_id}")
        
        # 获取GO Biological Process结果
        time.sleep(1)
        go_bp_url = f"{ENRICHR_URL}/export?userListId={user_list_id}&filename=GO_BP&backgroundType=GO_Biological_Process_2023"
        response = requests.get(go_bp_url, timeout=30)
        
        if response.status_code == 200:
            go_bp = response.json().get('GO_Biological_Process_2023', [])[:10]
            results['GO_BP'] = go_bp
            print(f"\nGO Biological Process (Top 5):")
            for item in go_bp[:5]:
                print(f"  {item[1]}: p={item[2]:.2e}, genes={item[5]}")
        
        # 获取KEGG结果
        time.sleep(1)
        kegg_url = f"{ENRICHR_URL}/export?userListId={user_list_id}&filename=KEGG&backgroundType=KEGG_2021_Human"
        response = requests.get(kegg_url, timeout=30)
        
        if response.status_code == 200:
            kegg = response.json().get('KEGG_2021_Human', [])[:10]
            results['KEGG'] = kegg
            print(f"\nKEGG Pathways (Top 5):")
            for item in kegg[:5]:
                print(f"  {item[1]}: p={item[2]:.2e}, genes={item[5]}")
        
        return results
        
    except Exception as e:
        print(f"Enrichr API error: {e}")
        return None


def hypergeometric_test(conn, G, acupoint_nodes, gene_nodes):
    """超几何分布检验：评估穴位-基因共现是否显著高于随机期望"""
    print("\n" + "=" * 60)
    print("4. 超几何分布显著性检验")
    print("=" * 60)
    
    cursor = conn.cursor()
    cursor.execute('SELECT COUNT(*) FROM articles')
    total_articles = cursor.fetchone()[0]
    
    cursor.execute('''
        SELECT g.gene_symbol, COUNT(DISTINCT ag.article_id) as cnt
        FROM genes g JOIN article_gene ag ON g.id = ag.gene_id
        WHERE g.quality_flag IN ('valid', 'mapped_to_official')
        GROUP BY g.id
    ''')
    gene_freq = {row[0]: row[1] for row in cursor.fetchall()}
    
    cursor.execute('''
        SELECT a.name_cn, COUNT(DISTINCT aa.article_id) as cnt
        FROM acupoints a JOIN article_acupoint aa ON a.id = aa.acupoint_id
        GROUP BY a.id
    ''')
    acupoint_freq = {row[0]: row[1] for row in cursor.fetchall()}
    
    # 对高频共现对进行超几何检验
    cursor.execute('''
        SELECT a.name_cn, g.gene_symbol, COUNT(*) as cnt
        FROM acupoints a
        JOIN article_acupoint aa ON a.id = aa.acupoint_id
        JOIN article_gene ag ON aa.article_id = ag.article_id
        JOIN genes g ON ag.gene_id = g.id
        WHERE g.quality_flag IN ('valid', 'mapped_to_official')
        GROUP BY a.id, g.id
        HAVING cnt >= 3
        ORDER BY cnt DESC
        LIMIT 20
    ''')
    
    results = []
    print(f"\n{'穴位':<10} {'基因':<12} {'共现':<6} {'期望':<8} {'p值':<12} {'显著性'}")
    print("-" * 70)
    
    for apt, gene, cooccur in cursor.fetchall():
        # 超几何分布参数
        N = total_articles  # 总体本量
        K = gene_freq.get(gene, 0)  # 基因出现的文献数
        n = acupoint_freq.get(apt, 0)  # 穴位出现的文献数
        k = cooccur  # 共现文献数
        
        if n == 0 or K == 0:
            continue
        
        # 期望共现数 (随机期望)
        expected = (K * n) / N
        
        # 超几何检验 p-value (共现数 >= k 的概率)
        # P(X >= k) = 1 - P(X < k) = 1 - cdf(k-1)
        try:
            p_value = 1 - stats.hypergeom.cdf(k - 1, N, K, n)
        except:
            p_value = 1.0
        
        significance = ""
        if p_value < 0.001:
            significance = "***"
        elif p_value < 0.01:
            significance = "**"
        elif p_value < 0.05:
            significance = "*"
        else:
            significance = "ns"
        
        results.append({
            'acupoint': apt,
            'gene': gene,
            'cooccur': cooccur,
            'expected': round(expected, 2),
            'p_value': p_value,
            'significance': significance
        })
        
        print(f"{apt:<10} {gene:<12} {cooccur:<6} {expected:<8.2f} {p_value:<12.2e} {significance}")
    
    return results


def visualize_network(G, partition, acupoint_nodes, gene_nodes, output_path):
    """可视化网络及社区"""
    fig, ax = plt.subplots(figsize=(16, 16))
    pos = nx.spring_layout(G, k=2.8, iterations=100, seed=42)
    
    # 根据社区分配颜色
    communities = defaultdict(list)
    for node, comm_id in partition.items():
        communities[comm_id].append(node)
    
    community_colors = plt.cm.Set3(np.linspace(0, 1, len(communities)))
    node_colors = {}
    for comm_id, nodes in communities.items():
        color = community_colors[comm_id % len(community_colors)]
        for node in nodes:
            node_colors[node] = color
    
    # 绘制边
    edges = G.edges(data=True)
    weights = [d['weight'] for (u, v, d) in edges]
    nx.draw_networkx_edges(G, pos, width=[max(0.3, w/10) for w in weights], alpha=0.25, edge_color='#888', ax=ax)
    
    # 绘制穴位节点（圆形，较大）
    apt_list = list(acupoint_nodes)
    nx.draw_networkx_nodes(G, pos, nodelist=apt_list, 
                           node_color=[node_colors.get(n, '#e74c3c') for n in apt_list],
                           node_size=900, edgecolors='white', linewidths=2, ax=ax)
    
    # 绘制基因节点（根据社区着色）
    gene_list = list(gene_nodes)
    nx.draw_networkx_nodes(G, pos, nodelist=gene_list,
                           node_color=[node_colors.get(n, '#3498db') for n in gene_list],
                           node_size=[150 + G.degree(n) * 60 for n in gene_list],
                           alpha=0.8, edgecolors='white', linewidths=1, ax=ax)
    
    # 标签
    nx.draw_networkx_labels(G, pos, font_size=8, ax=ax)
    
    ax.set_title('Acupoint-Gene Co-occurrence Network with Louvain Communities', 
                 fontsize=16, fontweight='bold')
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"\n网络社区图已保存: {output_path}")


def export_results(topology_results, cluster_results, enrichment_results, hypergeo_results):
    """导出所有分析结果"""
    # 网络拓扑
    with open(f'{OUTPUT_DIR}/network_topology.json', 'w', encoding='utf-8') as f:
        json.dump(topology_results, f, ensure_ascii=False, indent=2)
    
    # 聚类
    with open(f'{OUTPUT_DIR}/cluster_communities.json', 'w', encoding='utf-8') as f:
        json.dump(cluster_results, f, ensure_ascii=False, indent=2)
    
    # 富集分析
    if enrichment_results:
        with open(f'{OUTPUT_DIR}/enrichment_results.json', 'w', encoding='utf-8') as f:
            json.dump(enrichment_results, f, ensure_ascii=False, indent=2)
    
    # 超几何检验
    with open(f'{OUTPUT_DIR}/hypergeometric_test.csv', 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.DictWriter(f, fieldnames=['acupoint', 'gene', 'cooccur', 'expected', 'p_value', 'significance'])
        writer.writeheader()
        for row in hypergeo_results:
            writer.writerow(row)
    
    print(f"\n所有分析结果已导出至 {OUTPUT_DIR}/")


def main():
    print("=" * 60)
    print("针灸穴位-基因关联数据库 高级数据分析")
    print("=" * 60)
    
    conn = sqlite3.connect('output/acupoint_gene.db')
    
    # 1. 构建网络
    G, acupoint_nodes, gene_nodes, edges = build_network(conn)
    print(f"\n网络构建完成: {G.number_of_nodes()} 节点, {G.number_of_edges()} 边")
    
    # 2. 网络拓扑分析
    topology_results = network_topology_analysis(G, acupoint_nodes, gene_nodes)
    
    # 3. 聚类分析
    cluster_results, partition = cluster_analysis(G, acupoint_nodes, gene_nodes)
    
    # 4. 可视化网络社区
    visualize_network(G, partition, acupoint_nodes, gene_nodes, 
                     f'{OUTPUT_DIR}/fig_network_communities.png')
    
    # 5. 富集分析（对高频基因）
    top_genes = [g for g, _ in topology_results['degree_centrality']['top_genes'][:30]]
    enrichment_results = enrichment_analysis(top_genes, "Top genes from acupoint-gene network")
    
    # 6. 超几何显著性检验
    hypergeo_results = hypergeometric_test(conn, G, acupoint_nodes, gene_nodes)
    
    # 7. 导出结果
    export_results(topology_results, cluster_results, enrichment_results, hypergeo_results)
    
    conn.close()
    
    print("\n" + "=" * 60)
    print("分析完成！")
    print("=" * 60)


if __name__ == '__main__':
    main()
