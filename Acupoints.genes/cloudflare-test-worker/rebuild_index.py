#!/usr/bin/env python3
"""
Safely rebuild src/index.js with cleaned DB data and new analysis functions.
"""
import json

# Read original index.js
with open('src/index.js', 'r') as f:
    original_lines = f.readlines()

# Read cleaned DB data
with open('db_data.json', 'r') as f:
    db_data = json.load(f)

# Clean pmid_genes: remove gene families from pmid_genes 
# (they should only appear in gene_families data, not in per-article gene lists)
gene_families = {'AKT', 'PI3K', 'MAPK', 'WNT', 'JAK', 'ERK', 'MMP', 'BCL', 'P38', 'STAT'}
for pmid, genes in db_data.get('pmid_genes', {}).items():
    db_data['pmid_genes'][pmid] = [g for g in genes if g not in gene_families]

# Rebuild network_data from cleaned pmid_genes and pmid_acupoints
network_counts = {}
for pmid, apts in db_data.get('pmid_acupoints', {}).items():
    genes = db_data.get('pmid_genes', {}).get(pmid, [])
    for apt in apts:
        apt_name = apt.get('name', '')
        for gene in genes:
            key = (apt_name, gene)
            network_counts[key] = network_counts.get(key, 0) + 1

db_data['network_data'] = [
    {"acupoint": apt, "gene": gene, "count": count}
    for (apt, gene), count in sorted(network_counts.items(), key=lambda x: -x[1])
    if count >= 1
]

# Update stats
db_data['stats']['network_edges'] = len(db_data['network_data'])

# Normalize all PMID keys/values to strings for consistent JavaScript object lookup
db_data['acupoint_pmids'] = {k: [str(p) for p in v] for k, v in db_data.get('acupoint_pmids', {}).items()}
db_data['gene_pmids'] = {k: [str(p) for p in v] for k, v in db_data.get('gene_pmids', {}).items()}

# Convert DB data to a single line for inline embedding
db_json = json.dumps(db_data, ensure_ascii=False, separators=(',', ':'))

# New analysis functions to insert after getGeneFunction (after line 23, before getAcupointInfo)
new_functions = r'''
// ========================================================================
// 高级分析函数 (v3.0)
// ========================================================================
function getNetworkTopology() {
  const G = {};
  const degrees = {};
  const betweenness = {};
  const acupointNodes = new Set();
  const geneNodes = new Set();

  for (const edge of DB.network_data) {
    const { acupoint, gene, count } = edge;
    if (!G[acupoint]) G[acupoint] = {};
    if (!G[gene]) G[gene] = {};
    G[acupoint][gene] = count;
    G[gene][acupoint] = count;
    degrees[acupoint] = (degrees[acupoint] || 0) + 1;
    degrees[gene] = (degrees[gene] || 0) + 1;
    acupointNodes.add(acupoint);
    geneNodes.add(gene);
  }

  const totalNodes = Object.keys(G).length;

  // Degree centrality
  const degreeCentrality = {};
  for (const [node, deg] of Object.entries(degrees)) {
    degreeCentrality[node] = totalNodes > 1 ? deg / (totalNodes - 1) : 0;
  }

  // Betweenness centrality (approximate for bipartite)
  for (const node of Object.keys(G)) {
    betweenness[node] = 0;
  }
  for (const s of Object.keys(G)) {
    const queue = [s];
    const visited = new Set([s]);
    const parent = {};
    const dist = { [s]: 0 };
    let idx = 0;
    while (idx < queue.length) {
      const u = queue[idx++];
      for (const v of Object.keys(G[u] || {})) {
        if (!visited.has(v)) {
          visited.add(v);
          dist[v] = dist[u] + 1;
          parent[v] = [u];
          queue.push(v);
        } else if (dist[v] === dist[u] + 1) {
          parent[v].push(u);
        }
      }
    }
    const dependency = {};
    for (const n of queue) dependency[n] = 0;
    while (queue.length > 1) {
      const w = queue.pop();
      if (parent[w]) {
        for (const v of parent[w]) {
          dependency[v] += (1 + dependency[w]) / parent[w].length;
        }
      }
      if (w !== s) betweenness[w] += dependency[w];
    }
  }
  const maxBetweenness = Math.max(...Object.values(betweenness), 1);
  for (const n of Object.keys(betweenness)) {
    betweenness[n] = betweenness[n] / maxBetweenness;
  }

  const aptDegree = [...acupointNodes].map(n => ({ node: n, centrality: degreeCentrality[n] || 0 })).sort((a, b) => b.centrality - a.centrality);
  const geneDegree = [...geneNodes].map(n => ({ node: n, centrality: degreeCentrality[n] || 0 })).sort((a, b) => b.centrality - a.centrality);
  const geneBetweenness = [...geneNodes].map(n => ({ node: n, centrality: betweenness[n] || 0 })).sort((a, b) => b.centrality - a.centrality);

  return {
    network_stats: { total_nodes: totalNodes, total_edges: DB.network_data.length, density: totalNodes > 1 ? (DB.network_data.length * 2) / (totalNodes * (totalNodes - 1)) : 0 },
    degree_centrality: { acupoints: aptDegree.slice(0, 10), genes: geneDegree.slice(0, 15) },
    betweenness_centrality: { genes: geneBetweenness.slice(0, 15) }
  };
}

function getClusterCommunities() {
  // Louvain algorithm simplified
  const G = {};
  for (const edge of DB.network_data) {
    const { acupoint, gene, count } = edge;
    if (!G[acupoint]) G[acupoint] = {};
    if (!G[gene]) G[gene] = {};
    G[acupoint][gene] = count;
    G[gene][acupoint] = count;
  }

  const nodes = Object.keys(G);
  const partition = {};
  for (let i = 0; i < nodes.length; i++) partition[nodes[i]] = i;

  function getModularity(part) {
    let Q = 0;
    const m = DB.network_data.reduce((s, e) => s + e.count, 0);
    const commTotal = {};
    const commInternal = {};
    for (const edge of DB.network_data) {
      const c1 = part[edge.acupoint];
      const c2 = part[edge.gene];
      commTotal[c1] = (commTotal[c1] || 0) + edge.count;
      commTotal[c2] = (commTotal[c2] || 0) + edge.count;
      if (c1 === c2) {
        commInternal[c1] = (commInternal[c1] || 0) + edge.count * 2;
      }
    }
    for (const c of Object.keys(commInternal)) {
      Q += (commInternal[c] / (2 * m)) - Math.pow(commTotal[c] / (2 * m), 2);
    }
    return Q;
  }

  // Greedy optimization (simplified)
  let improved = true;
  let iterations = 0;
  while (improved && iterations < 10) {
    improved = false;
    iterations++;
    for (const node of nodes) {
      const currentComm = partition[node];
      const neighborComms = {};
      for (const neighbor of Object.keys(G[node] || {})) {
        const comm = partition[neighbor];
        neighborComms[comm] = (neighborComms[comm] || 0) + G[node][neighbor];
      }
      let bestComm = currentComm;
      let bestGain = 0;
      for (const [comm, weight] of Object.entries(neighborComms)) {
        if (comm !== currentComm) {
          const gain = weight;
          if (gain > bestGain) {
            bestGain = gain;
            bestComm = parseInt(comm);
          }
        }
      }
      if (bestComm !== currentComm) {
        partition[node] = bestComm;
        improved = true;
      }
    }
  }

  // Renumber communities
  const commMap = {};
  let nextId = 0;
  const finalPartition = {};
  for (const node of nodes) {
    const oldId = partition[node];
    if (!(oldId in commMap)) commMap[oldId] = nextId++;
    finalPartition[node] = commMap[oldId];
  }

  const communities = {};
  for (const [node, commId] of Object.entries(finalPartition)) {
    if (!communities[commId]) communities[commId] = { acupoints: [], genes: [] };
    if (DB.acupoint_info[node]) communities[commId].acupoints.push(node);
    else communities[commId].genes.push(node);
  }

  return { num_communities: nextId, communities: Object.entries(communities).map(([id, c]) => ({ community_id: parseInt(id), acupoints: c.acupoints, genes: c.genes, size: c.acupoints.length + c.genes.length })).sort((a, b) => b.size - a.size) };
}

function getHypergeometricTest() {
  const totalArticles = DB.stats.articles;
  const geneFreq = {};
  const acupointFreq = {};

  for (const [pmid, genes] of Object.entries(DB.pmid_genes)) {
    for (const g of genes) {
      geneFreq[g] = (geneFreq[g] || 0) + 1;
    }
  }
  for (const [pmid, apts] of Object.entries(DB.pmid_acupoints)) {
    for (const a of apts) {
      acupointFreq[a.name] = (acupointFreq[a.name] || 0) + 1;
    }
  }

  const results = [];
  for (const edge of DB.network_data.slice(0, 20)) {
    const { acupoint, gene, count } = edge;
    const K = geneFreq[gene] || 0;
    const n = acupointFreq[acupoint] || 0;
    const expected = (K * n) / totalArticles;

    // Approximate p-value using binomial approximation
    const p = K / totalArticles;
    const pValue = 1 - binomialCDF(count - 1, n, p);

    let sig = "ns";
    if (pValue < 0.001) sig = "***";
    else if (pValue < 0.01) sig = "**";
    else if (pValue < 0.05) sig = "*";

    results.push({ acupoint, gene, cooccur: count, expected: Math.round(expected * 100) / 100, p_value: pValue, significance: sig });
  }
  return results;
}

function binomialCDF(k, n, p) {
  let cdf = 0;
  for (let i = 0; i <= k; i++) {
    cdf += binomialPMF(i, n, p);
  }
  return cdf;
}

function binomialPMF(k, n, p) {
  if (k < 0 || k > n) return 0;
  const coeff = factorial(n) / (factorial(k) * factorial(n - k));
  return coeff * Math.pow(p, k) * Math.pow(1 - p, n - k);
}

function factorial(n) {
  if (n <= 1) return 1;
  let result = 1;
  for (let i = 2; i <= n; i++) result *= i;
  return result;
}

'''

# New render functions to insert after renderAbout (before export default)
render_topology = r'''
function renderTopology() {
  const data = getNetworkTopology();
  const stats = data.network_stats;
  const degreeApt = data.degree_centrality.acupoints;
  const degreeGene = data.degree_centrality.genes;
  const betweennessGene = data.betweenness_centrality.genes;

  return COMMON_HEAD("网络拓扑分析") + NAV + `
<div class="fade-in">
  <h1 style="margin-bottom:8px;font-size:1.5em">🕸️ 网络拓扑分析</h1>
  <p style="color:var(--text-light);margin-bottom:20px">基于穴位-基因共现网络的拓扑结构分析</p>
  <div class="stat-cards">
    <div class="stat-card" style="border-top-color:#1a5fb4"><div class="number">${stats.total_nodes}</div><div class="label">网络节点</div></div>
    <div class="stat-card" style="border-top-color:#26a269"><div class="number">${stats.total_edges}</div><div class="label">网络边数</div></div>
    <div class="stat-card" style="border-top-color:#c061cb"><div class="number">${Math.round(stats.density * 10000) / 10000}</div><div class="label">网络密度</div></div>
  </div>
  <div style="display:grid;grid-template-columns:1fr 1fr;gap:20px;margin-top:20px">
    <div class="card"><div class="card-title">📍 穴位度中心性 Top 10</div>
      <table class="data-table"><thead><tr><th>穴位</th><th>度中心性</th></tr></thead><tbody>
        ${degreeApt.map(d => '<tr><td>' + d.node + '</td><td>' + (Math.round(d.centrality * 10000) / 10000) + '</td></tr>').join('')}
      </tbody></table>
    </div>
    <div class="card"><div class="card-title">🧬 基因度中心性 Top 15</div>
      <table class="data-table"><thead><tr><th>基因</th><th>度中心性</th></tr></thead><tbody>
        ${degreeGene.map(d => '<tr><td>' + d.node + '</td><td>' + (Math.round(d.centrality * 10000) / 10000) + '</td></tr>').join('')}
      </tbody></table>
    </div>
  </div>
  <div class="card" style="margin-top:20px"><div class="card-title">🌉 基因中介中心性 Top 15</div>
    <table class="data-table"><thead><tr><th>基因</th><th>中介中心性</th></tr></thead><tbody>
      ${betweennessGene.map(d => '<tr><td>' + d.node + '</td><td>' + (Math.round(d.centrality * 10000) / 10000) + '</td></tr>').join('')}
    </tbody></table>
  </div>
</div>
` + FOOTER;
}

'''

render_clusters = r'''
function renderClusters() {
  const data = getClusterCommunities();
  const communities = data.communities;

  return COMMON_HEAD("聚类分析") + NAV + `
<div class="fade-in">
  <h1 style="margin-bottom:8px;font-size:1.5em">🔍 Louvain 聚类分析</h1>
  <p style="color:var(--text-light);margin-bottom:20px">基于模块度优化的社区发现算法，识别穴位-基因功能模块</p>
  <div class="stat-cards">
    <div class="stat-card" style="border-top-color:#1a5fb4"><div class="number">${data.num_communities}</div><div class="label">社区数量</div></div>
    <div class="stat-card" style="border-top-color:#26a269"><div class="number">${communities.length}</div><div class="label">有效社区</div></div>
  </div>
  <div style="margin-top:20px">
    ${communities.map((c, idx) => '\n      <div class="card" style="margin-bottom:16px">\n        <div class="card-title">社区 ' + c.community_id + ' (大小: ' + c.size + ')</div>\n        <p style="margin:8px 0"><strong>穴位:</strong> ' + (c.acupoints.join('、') || '无') + '</p>\n        <p style="margin:8px 0"><strong>基因:</strong> ' + c.genes.slice(0, 30).join('、') + (c.genes.length > 30 ? ' ...' : '') + '</p>\n      </div>\n    ').join('')}
  </div>
</div>
` + FOOTER;
}

'''

# Build the new file
new_lines = []

# Lines 1-14 (up to and including the DB data comment)
new_lines.extend(original_lines[:14])

# Line 15: new DB data
new_lines.append('const DB = ' + db_json + ';\n')

# Line 16+: rest of original up to line 23 (getGeneFunction end)
new_lines.extend(original_lines[15:23])

# Insert new analysis functions after getGeneFunction (after line 23)
new_lines.append(new_functions)

# Lines 24+ of original (from getAcupointInfo onward)
new_lines.extend(original_lines[23:])

# Now find the position of renderAbout end and insert new render functions before export default
# We need to re-join and re-split to properly track line numbers
content = ''.join(new_lines)
lines_after_insert = content.split('\n')
lines_after_insert = [line + '\n' for line in lines_after_insert[:-1]] + ([lines_after_insert[-1]] if lines_after_insert[-1] else [])

# Find where renderAbout ends
in_func = False
brace_count = 0
render_about_end = -1
for i, line in enumerate(lines_after_insert):
    if 'function renderAbout' in line:
        in_func = True
        brace_count = 0
    if in_func:
        for char in line:
            if char == '{':
                brace_count += 1
            elif char == '}':
                brace_count -= 1
                if brace_count == 0:
                    render_about_end = i
                    in_func = False
                    break

print(f"renderAbout ends at line {render_about_end + 1}")

# Find export default line
export_default_line = -1
for i, line in enumerate(lines_after_insert):
    if line.strip() == 'export default {':
        export_default_line = i
        break

print(f"export default at line {export_default_line + 1}")

# Insert render functions after renderAbout, before export default
final_lines = []
final_lines.extend(lines_after_insert[:render_about_end + 1])
final_lines.append('\n')
final_lines.append(render_topology)
final_lines.append(render_clusters)
final_lines.extend(lines_after_insert[render_about_end + 1:])

# Update NAV to add topology and clusters links
content2 = ''.join(final_lines)
content2 = content2.replace(
    '    <a href="/network" id="nav-network">共现网络</a>\n    <a href="/about" id="nav-about">算法说明</a>',
    '    <a href="/network" id="nav-network">共现网络</a>\n    <a href="/topology" id="nav-topology">网络拓扑</a>\n    <a href="/clusters" id="nav-clusters">聚类分析</a>\n    <a href="/about" id="nav-about">算法说明</a>'
)
content2 = content2.replace(
    '    <a href="/network" id="side-network">🕸️ 共现网络</a>\n    <a href="/about" id="side-about">ℹ️ 算法说明</a>',
    '    <a href="/network" id="side-network">🕸️ 共现网络</a>\n    <a href="/topology" id="side-topology">🌐 网络拓扑</a>\n    <a href="/clusters" id="side-clusters">🔍 聚类分析</a>\n    <a href="/about" id="side-about">ℹ️ 算法说明</a>'
)

# Add new API routes in fetch handler
api_insertion = '''      if (path === "/api/topology") {
        return jsonResponse(getNetworkTopology());
      }
      if (path === "/api/clusters") {
        return jsonResponse(getClusterCommunities());
      }
      if (path === "/api/hypergeo") {
        return jsonResponse(getHypergeometricTest());
      }
'''

content2 = content2.replace(
    '      if (path === "/api/network") {',
    api_insertion + '      if (path === "/api/network") {'
)

# Add new page routes in fetch handler
page_insertion = '''      if (path === "/topology") {
        return htmlResponse(renderTopology());
      }
      if (path === "/clusters") {
        return htmlResponse(renderClusters());
      }
'''

content2 = content2.replace(
    '      if (path === "/about") {',
    page_insertion + '      if (path === "/about") {'
)

with open('src/index.js', 'w') as f:
    f.write(content2)

print(f"Rebuilt src/index.js ({len(content2)} bytes)")

# Verify syntax with node (if available)
import subprocess
try:
    result = subprocess.run(['node', '--check', 'src/index.js'], capture_output=True, text=True, timeout=10)
    if result.returncode == 0:
        print("✓ Syntax check passed")
    else:
        print(f"✗ Syntax error:\n{result.stderr}")
except FileNotFoundError:
    print("node not available for syntax check")
except Exception as e:
    print(f"Syntax check error: {e}")
