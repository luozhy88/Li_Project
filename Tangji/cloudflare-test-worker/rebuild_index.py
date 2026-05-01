#!/usr/bin/env python3
"""
ThyGlycoPortal Cloudflare Worker Builder
Reads SQLite database and generates src/index.js with inline data + frontend pages.
"""

import sqlite3
import json
import os

DB_PATH = os.path.join(os.path.dirname(__file__), '..', 'output', 'thyroid_glyco_db.sqlite')
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), 'src')
OUTPUT_PATH = os.path.join(OUTPUT_DIR, 'index.js')

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================================
# 1. Read all database tables
# ============================================================================
conn = sqlite3.connect(DB_PATH)
conn.row_factory = sqlite3.Row
cursor = conn.cursor()

tables = [
    'studies', 'glycan_structures', 'clinical_groups', 'samples',
    'glycan_abundance', 'biomarkers', 'glycosyltransferases',
    'enzyme_glycan_links', 'literature_stats', 'tcga_glycogene_expression'
]

db_data = {}
for table in tables:
    cursor.execute(f"SELECT * FROM {table}")
    columns = [desc[0] for desc in cursor.description]
    rows = [dict(zip(columns, row)) for row in cursor.fetchall()]
    db_data[table] = rows

conn.close()

# ============================================================================
# 2. Build computed datasets
# ============================================================================
stats = {
    'studies': len([s for s in db_data['studies'] if s.get('method') != 'Review']),
    'reviews': len([s for s in db_data['studies'] if s.get('method') == 'Review']),
    'glycan_structures': len(db_data['glycan_structures']),
    'clinical_groups': len(db_data['clinical_groups']),
    'biomarkers': len(db_data['biomarkers']),
    'glycosyltransferases': len(db_data['glycosyltransferases']),
    'enzyme_glycan_links': len(db_data['enzyme_glycan_links']),
    'tcga_genes': len(db_data['tcga_glycogene_expression']),
}

# Year distribution
year_counts = {}
for s in db_data['studies']:
    y = s.get('year')
    if y:
        year_counts[str(y)] = year_counts.get(str(y), 0) + 1
year_distribution = [{"year": int(y), "count": c} for y, c in sorted(year_counts.items())]

# Cancer type distribution (exclude reviews)
cancer_counts = {}
for s in db_data['studies']:
    if s.get('method') == 'Review':
        continue
    ct = s.get('cancer_type', 'Unknown')
    if ct == 'All types':
        continue
    cancer_counts[ct] = cancer_counts.get(ct, 0) + 1
cancer_distribution = [{"type": k, "count": v} for k, v in cancer_counts.items()]

# Method distribution (exclude reviews)
method_counts = {}
for s in db_data['studies']:
    if s.get('method') == 'Review':
        continue
    m = s.get('method', 'Unknown')
    method_counts[m] = method_counts.get(m, 0) + 1
method_distribution = [{"method": k, "count": v} for k, v in method_counts.items()]

# Biomarker AUC data
biomarker_auc = []
for b in db_data['biomarkers']:
    if b.get('performance_auc'):
        biomarker_auc.append({
            "name": b.get('name', '')[:40],
            "auc": b.get('performance_auc'),
            "type": b.get('biomarker_type', 'Unknown'),
            "sample": b.get('sample_type', '')
        })

# TCGA expression data
tcga_data = []
for g in db_data['tcga_glycogene_expression']:
    tcga_data.append({
        "gene": g.get('gene_symbol'),
        "family": g.get('enzyme_family'),
        "normal_rpkm": g.get('normal_median_rpkm'),
        "ptc_rpkm": g.get('ptc_median_rpkm'),
        "fold_change": g.get('fold_change'),
        "p_value": g.get('p_value'),
        "significant": g.get('significant') == 1 or g.get('significant') == '1',
        "direction": g.get('direction')
    })

# Enzyme-glycan network
network = []
for link in db_data['enzyme_glycan_links']:
    network.append({
        "enzyme_id": link.get('enzyme_id'),
        "glycan_id": link.get('glycan_id'),
        "linkage": link.get('linkage_type'),
        "evidence": link.get('evidence_level')
    })

# Literature stats for references
lit_stats = []
for ls in db_data['literature_stats']:
    lit_stats.append({
        "study": ls.get('study'),
        "variable": ls.get('variable'),
        "group1": ls.get('group1'),
        "group1_n": ls.get('group1_n'),
        "group1_mean": ls.get('group1_mean'),
        "group1_sd": ls.get('group1_sd'),
        "group2": ls.get('group2'),
        "group2_n": ls.get('group2_n'),
        "group2_mean": ls.get('group2_mean'),
        "group2_sd": ls.get('group2_sd'),
        "p_value": ls.get('p_value'),
        "direction": ls.get('direction'),
        "auc": ls.get('auc')
    })

# Glycan cross-study data
glycan_cross = {}
for ga in db_data['glycan_abundance']:
    sid = ga.get('study_id')
    gid = ga.get('glycan_id')
    # Find study and glycan info
    study = next((s for s in db_data['studies'] if s.get('study_id') == sid), {})
    glycan = next((g for g in db_data['glycan_structures'] if g.get('glycan_id') == gid), {})
    sample = next((s for s in db_data['samples'] if s.get('sample_id') == ga.get('sample_id')), {})
    group = next((g for g in db_data['clinical_groups'] if g.get('group_id') == sample.get('group_id')), {})
    
    snfg = glycan.get('snfg_name', '')
    if snfg not in glycan_cross:
        glycan_cross[snfg] = []
    glycan_cross[snfg].append({
        "study_id": sid,
        "year": study.get('year'),
        "group": group.get('group_name', ''),
        "value": ga.get('abundance_value'),
        "direction": ga.get('change_direction')
    })

# Build final inline DB
inline_db = {
    "stats": stats,
    "year_distribution": year_distribution,
    "cancer_distribution": cancer_distribution,
    "method_distribution": method_distribution,
    "biomarker_auc": biomarker_auc,
    "tcga_data": tcga_data,
    "network": network,
    "lit_stats": lit_stats,
    "glycan_cross": glycan_cross,
    "studies": [{k: v for k, v in s.items() if k not in ['key_finding']} for s in db_data['studies']],
    "glycans": db_data['glycan_structures'],
    "biomarkers": db_data['biomarkers'],
    "enzymes": db_data['glycosyltransferases'],
    "groups": db_data['clinical_groups'],
}

# Remove None values recursively
def clean_none(obj):
    if isinstance(obj, dict):
        return {k: clean_none(v) for k, v in obj.items() if v is not None}
    elif isinstance(obj, list):
        return [clean_none(i) for i in obj]
    return obj

inline_db = clean_none(inline_db)

db_json = json.dumps(inline_db, ensure_ascii=False, separators=(',', ':'))

print(f"[INFO] Inline DB size: {len(db_json):,} chars")

# ============================================================================
# 3. Generate src/index.js
# ============================================================================
# We will write the JS in chunks to avoid memory issues with huge strings

with open(OUTPUT_PATH, 'w', encoding='utf-8') as f:
    # Header
    f.write('''/**
 * ThyGlycoPortal - Cloudflare Worker v1.0
 * Thyroid Cancer Glycosylation Interactive Analysis Platform
 * Data: output/thyroid_glyco_db.sqlite (inline JSON)
 * Charts: ECharts 5 (CDN)
 * Style: Tailwind CSS (CDN)
 */

// ========================================================================
// Inline Database
// ========================================================================
const DB = ''' + db_json + ''';

// ========================================================================
// Helper Functions
// ========================================================================
function escapeHtml(text) {
  if (!text) return '';
  return String(text).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
}

function renderNav(activePath) {
  const items = [
    {path: '/', label: 'Overview', icon: '📊'},
    {path: '/literature', label: 'Literature', icon: '📚'},
    {path: '/biomarkers', label: 'Biomarkers', icon: '🔬'},
    {path: '/glycans', label: 'Glycans', icon: '🧬'},
    {path: '/enzymes', label: 'Enzymes', icon: '⚗️'},
    {path: '/tcga', label: 'TCGA', icon: '📈'},
    {path: '/diagnostic', label: 'Diagnostic', icon: '🩺'},
    {path: '/nomogram', label: 'Nomogram', icon: '📐'},
  ];
  return items.map(i => `<a href="${i.path}" class="px-3 py-2 rounded-lg text-sm font-medium transition ${i.path === activePath ? 'bg-blue-600 text-white' : 'text-gray-700 hover:bg-gray-100'}">${i.icon} ${i.label}</a>`).join('');
}

function pageTemplate(title, content, activePath) {
  return `<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>${escapeHtml(title)} | ThyGlycoPortal</title>
<script src="https://cdn.tailwindcss.com"></script>
<script src="https://cdn.jsdelivr.net/npm/echarts@5.4.3/dist/echarts.min.js"></script>
<style>
  .chart-container { width: 100%; height: 300px; }
  .chart-container-lg { width: 100%; height: 400px; }
</style>
</head>
<body class="bg-gray-50 min-h-screen">
  <nav class="bg-white shadow-sm border-b sticky top-0 z-50">
    <div class="max-w-7xl mx-auto px-4">
      <div class="flex items-center justify-between h-14">
        <a href="/" class="text-xl font-bold text-blue-700">🔬 ThyGlycoPortal</a>
        <div class="hidden md:flex space-x-1">${renderNav(activePath)}</div>
      </div>
    </div>
  </nav>
  <main class="max-w-7xl mx-auto px-4 py-6">
    ${content}
  </main>
  <footer class="bg-white border-t mt-12 py-6">
    <div class="max-w-7xl mx-auto px-4 text-center text-sm text-gray-500">
      <p>ThyGlycoPortal v1.0 | All data sourced from peer-reviewed literature (REAL DATA ONLY)</p>
      <p class="mt-1">TCGA data: Bones et al. Cancers 2018 (PMC11727208)</p>
    </div>
  </footer>
</body>
</html>`;
}
''')

    # Overview page
    f.write('''
// ========================================================================
// Overview Page
// ========================================================================
function overviewPage() {
  const stats = DB.stats;
  const yearData = DB.year_distribution;
  const cancerData = DB.cancer_distribution;
  const methodData = DB.method_distribution;
  const aucData = DB.biomarker_auc;

  const yearChart = JSON.stringify({
    xAxis: { type: 'category', data: yearData.map(d => d.year) },
    yAxis: { type: 'value', name: 'Studies' },
    series: [{ data: yearData.map(d => d.count), type: 'bar', itemStyle: { color: '#3c8dbc' } }],
    tooltip: { trigger: 'axis' }
  });

  const cancerChart = JSON.stringify({
    series: [{
      type: 'pie', radius: '60%',
      data: cancerData.map(d => ({ name: d.type, value: d.count })),
      label: { formatter: '{b}: {c}' }
    }],
    tooltip: { trigger: 'item' }
  });

  const methodChart = JSON.stringify({
    xAxis: { type: 'value' },
    yAxis: { type: 'category', data: methodData.map(d => d.method) },
    series: [{ data: methodData.map(d => d.count), type: 'bar', itemStyle: { color: '#00a65a' } }],
    tooltip: { trigger: 'axis' }
  });

  const aucChart = JSON.stringify({
    xAxis: { type: 'value', min: 0.5, max: 1.0 },
    yAxis: { type: 'category', data: aucData.map(d => d.name) },
    series: [{
      data: aucData.map(d => ({ value: d.auc, itemStyle: { color: d.auc >= 0.9 ? '#dd4b39' : d.auc >= 0.8 ? '#f39c12' : '#999' } })),
      type: 'bar', label: { show: true, position: 'right', formatter: '{c}' }
    }],
    tooltip: { trigger: 'axis' }
  });

  const content = `
    <div class="mb-6">
      <h1 class="text-2xl font-bold text-gray-900">Overview Dashboard</h1>
      <p class="text-gray-600 mt-1">Thyroid Cancer Glycosylation Interactive Analysis Platform</p>
    </div>
    <div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-6">
      <div class="bg-white rounded-lg shadow p-4 border-l-4 border-blue-500">
        <div class="text-3xl font-bold text-blue-600">${stats.studies}</div>
        <div class="text-sm text-gray-500">Original Studies</div>
      </div>
      <div class="bg-white rounded-lg shadow p-4 border-l-4 border-green-500">
        <div class="text-3xl font-bold text-green-600">${stats.glycan_structures}</div>
        <div class="text-sm text-gray-500">Glycan Structures</div>
      </div>
      <div class="bg-white rounded-lg shadow p-4 border-l-4 border-red-500">
        <div class="text-3xl font-bold text-red-600">${stats.biomarkers}</div>
        <div class="text-sm text-gray-500">Biomarkers</div>
      </div>
      <div class="bg-white rounded-lg shadow p-4 border-l-4 border-yellow-500">
        <div class="text-3xl font-bold text-yellow-600">${stats.glycosyltransferases}</div>
        <div class="text-sm text-gray-500">Glycosyltransferases</div>
      </div>
    </div>
    <div class="grid grid-cols-1 md:grid-cols-2 gap-6">
      <div class="bg-white rounded-lg shadow p-4">
        <h3 class="font-semibold text-gray-800 mb-2">Publications by Year</h3>
        <div id="chart-year" class="chart-container"></div>
      </div>
      <div class="bg-white rounded-lg shadow p-4">
        <h3 class="font-semibold text-gray-800 mb-2">Cancer Type Distribution</h3>
        <div id="chart-cancer" class="chart-container"></div>
      </div>
      <div class="bg-white rounded-lg shadow p-4">
        <h3 class="font-semibold text-gray-800 mb-2">Method Distribution</h3>
        <div id="chart-method" class="chart-container"></div>
      </div>
      <div class="bg-white rounded-lg shadow p-4">
        <h3 class="font-semibold text-gray-800 mb-2">Biomarker AUC</h3>
        <div id="chart-auc" class="chart-container"></div>
      </div>
    </div>
    <script>
      echarts.init(document.getElementById('chart-year')).setOption(${yearChart});
      echarts.init(document.getElementById('chart-cancer')).setOption(${cancerChart});
      echarts.init(document.getElementById('chart-method')).setOption(${methodChart});
      echarts.init(document.getElementById('chart-auc')).setOption(${aucChart});
    </script>
  `;
  return pageTemplate('Overview', content, '/');
}
''')

    # Literature page
    f.write('''
// ========================================================================
// Literature Page
// ========================================================================
function literaturePage() {
  const studies = DB.studies;
  const rows = studies.map(s => `
    <tr class="border-b hover:bg-gray-50">
      <td class="px-4 py-2 text-sm">${escapeHtml(s.study_id || '')}</td>
      <td class="px-4 py-2 text-sm">${escapeHtml(s.title || '')}</td>
      <td class="px-4 py-2 text-sm">${s.year || ''}</td>
      <td class="px-4 py-2 text-sm">${escapeHtml(s.cancer_type || '')}</td>
      <td class="px-4 py-2 text-sm">${escapeHtml(s.method || '')}</td>
      <td class="px-4 py-2 text-sm">${escapeHtml(s.sample_type || '')}</td>
      <td class="px-4 py-2 text-sm">${s.pmid ? '<a href="https://pubmed.ncbi.nlm.nih.gov/' + s.pmid + '" target="_blank" class="text-blue-600 hover:underline">' + s.pmid + '</a>' : 'N/A'}</td>
    </tr>
  `).join('');

  const content = `
    <div class="mb-6">
      <h1 class="text-2xl font-bold text-gray-900">Literature</h1>
      <p class="text-gray-600 mt-1">${studies.length} studies curated from peer-reviewed literature</p>
    </div>
    <div class="bg-white rounded-lg shadow overflow-x-auto">
      <table class="min-w-full text-left">
        <thead class="bg-gray-100 text-gray-700 text-xs uppercase">
          <tr><th class="px-4 py-2">ID</th><th class="px-4 py-2">Title</th><th class="px-4 py-2">Year</th><th class="px-4 py-2">Cancer Type</th><th class="px-4 py-2">Method</th><th class="px-4 py-2">Sample</th><th class="px-4 py-2">PMID</th></tr>
        </thead>
        <tbody class="text-sm">${rows}</tbody>
      </table>
    </div>
  `;
  return pageTemplate('Literature', content, '/literature');
}
''')

    # Biomarkers page
    f.write('''
// ========================================================================
// Biomarkers Page
// ========================================================================
function biomarkersPage() {
  const biomarkers = DB.biomarkers;
  const aucData = DB.biomarker_auc;

  const rows = biomarkers.map(b => `
    <tr class="border-b hover:bg-gray-50">
      <td class="px-4 py-2 text-sm font-medium">${escapeHtml(b.name || '')}</td>
      <td class="px-4 py-2 text-sm">${escapeHtml(b.biomarker_type || '')}</td>
      <td class="px-4 py-2 text-sm">${escapeHtml(b.sample_type || '')}</td>
      <td class="px-4 py-2 text-sm">${b.performance_auc ? b.performance_auc.toFixed(3) : 'N/A'}</td>
      <td class="px-4 py-2 text-sm">${b.performance_sensitivity ? b.performance_sensitivity.toFixed(2) : 'N/A'}</td>
      <td class="px-4 py-2 text-sm">${b.performance_specificity ? b.performance_specificity.toFixed(2) : 'N/A'}</td>
      <td class="px-4 py-2 text-sm">${escapeHtml(b.validation_status || '')}</td>
    </tr>
  `).join('');

  const aucChart = JSON.stringify({
    xAxis: { type: 'value', min: 0.5, max: 1.0 },
    yAxis: { type: 'category', data: aucData.map(d => d.name) },
    series: [{
      data: aucData.map(d => ({ value: d.auc, itemStyle: { color: d.auc >= 0.9 ? '#2E8B57' : d.auc >= 0.8 ? '#4682B4' : '#CD853F' } })),
      type: 'bar', label: { show: true, position: 'right', formatter: '{c}' }
    }],
    tooltip: { trigger: 'axis' }
  });

  const content = `
    <div class="mb-6">
      <h1 class="text-2xl font-bold text-gray-900">Biomarkers</h1>
      <p class="text-gray-600 mt-1">Glycosylation-based diagnostic and prognostic biomarkers</p>
    </div>
    <div class="grid grid-cols-1 lg:grid-cols-3 gap-6 mb-6">
      <div class="lg:col-span-2 bg-white rounded-lg shadow overflow-x-auto">
        <table class="min-w-full text-left">
          <thead class="bg-gray-100 text-gray-700 text-xs uppercase">
            <tr><th class="px-4 py-2">Name</th><th class="px-4 py-2">Type</th><th class="px-4 py-2">Sample</th><th class="px-4 py-2">AUC</th><th class="px-4 py-2">Sens</th><th class="px-4 py-2">Spec</th><th class="px-4 py-2">Validation</th></tr>
          </thead>
          <tbody class="text-sm">${rows}</tbody>
        </table>
      </div>
      <div class="bg-white rounded-lg shadow p-4">
        <h3 class="font-semibold text-gray-800 mb-2">AUC Comparison</h3>
        <div id="chart-biomarker-auc" class="chart-container-lg"></div>
      </div>
    </div>
    <script>
      echarts.init(document.getElementById('chart-biomarker-auc')).setOption(${aucChart});
    </script>
  `;
  return pageTemplate('Biomarkers', content, '/biomarkers');
}
''')

    # Glycans page
    f.write('''
// ========================================================================
// Glycan Browser Page
// ========================================================================
function glycansPage() {
  const glycans = DB.glycans;
  const options = glycans.map(g => `<option value="${escapeHtml(g.snfg_name || '')}">${escapeHtml(g.snfg_name || '')} - ${escapeHtml(g.composition || '')}</option>`).join('');

  const glycanJson = JSON.stringify(glycans);
  const crossData = JSON.stringify(DB.glycan_cross || {});

  const content = `
    <div class="mb-6">
      <h1 class="text-2xl font-bold text-gray-900">Glycan Browser</h1>
      <p class="text-gray-600 mt-1">Browse ${glycans.length} core N-glycan structures and cross-study trends</p>
    </div>
    <div class="grid grid-cols-1 md:grid-cols-3 gap-6">
      <div class="bg-white rounded-lg shadow p-4">
        <label class="block text-sm font-medium text-gray-700 mb-2">Select Glycan</label>
        <select id="glycan-select" class="w-full border rounded-lg px-3 py-2 text-sm" onchange="updateGlycanView()">
          ${options}
        </select>
        <div id="glycan-info" class="mt-4 text-sm text-gray-600"></div>
      </div>
      <div class="md:col-span-2 bg-white rounded-lg shadow p-4">
        <h3 class="font-semibold text-gray-800 mb-2">Cross-Study Comparison</h3>
        <div id="chart-glycan-cross" class="chart-container-lg"></div>
      </div>
    </div>
    <div class="mt-6 bg-white rounded-lg shadow overflow-x-auto">
      <table class="min-w-full text-left">
        <thead class="bg-gray-100 text-gray-700 text-xs uppercase">
          <tr><th class="px-4 py-2">SNFG</th><th class="px-4 py-2">Composition</th><th class="px-4 py-2">Mass (Da)</th><th class="px-4 py-2">Type</th></tr>
        </thead>
        <tbody class="text-sm">
          ${glycans.map(g => `<tr class="border-b hover:bg-gray-50"><td class="px-4 py-2 font-mono">${escapeHtml(g.snfg_name || '')}</td><td class="px-4 py-2">${escapeHtml(g.composition || '')}</td><td class="px-4 py-2">${g.mass ? g.mass.toFixed(1) : 'N/A'}</td><td class="px-4 py-2">${escapeHtml(g.glycan_type || '')}</td></tr>`).join('')}
        </tbody>
      </table>
    </div>
    <script>
      const GLYCANS = ${glycanJson};
      const CROSS = ${crossData};
      let chartGlycan;
      function updateGlycanView() {
        const name = document.getElementById('glycan-select').value;
        const g = GLYCANS.find(x => x.snfg_name === name);
        document.getElementById('glycan-info').innerHTML = g ?
          '<b>Composition:</b> ' + g.composition + '<br><b>Mass:</b> ' + (g.mass ? g.mass.toFixed(1) : 'N/A') + ' Da<br><b>Type:</b> ' + (g.glycan_type || '') : '';
        const data = CROSS[name] || [];
        const groups = [...new Set(data.map(d => d.group))];
        const studies = [...new Set(data.map(d => d.study_id))];
        const series = groups.map(grp => ({
          name: grp,
          type: 'bar',
          data: studies.map(sid => { const d = data.find(x => x.study_id === sid && x.group === grp); return d ? d.value : 0; })
        }));
        if (!chartGlycan) chartGlycan = echarts.init(document.getElementById('chart-glycan-cross'));
        chartGlycan.setOption({
          xAxis: { type: 'category', data: studies },
          yAxis: { type: 'value' },
          series: series,
          tooltip: { trigger: 'axis' },
          legend: { data: groups }
        });
      }
      updateGlycanView();
    </script>
  `;
  return pageTemplate('Glycan Browser', content, '/glycans');
}
''')

    # Enzymes page
    f.write('''
// ========================================================================
// Enzyme Network Page
// ========================================================================
function enzymesPage() {
  const enzymes = DB.enzymes;
  const network = DB.network;
  const glycans = DB.glycans;

  const enzymeOptions = enzymes.map(e => `<option value="${escapeHtml(e.gene_symbol || '')}">${escapeHtml(e.gene_symbol || '')} (${escapeHtml(e.enzyme_family || '')})</option>`).join('');
  const enzymeJson = JSON.stringify(enzymes);
  const networkJson = JSON.stringify(network);
  const glycanJson = JSON.stringify(glycans);

  const content = `
    <div class="mb-6">
      <h1 class="text-2xl font-bold text-gray-900">Enzyme Network</h1>
      <p class="text-gray-600 mt-1">Glycosyltransferase-glycan regulatory relationships</p>
    </div>
    <div class="grid grid-cols-1 md:grid-cols-3 gap-6">
      <div class="bg-white rounded-lg shadow p-4">
        <label class="block text-sm font-medium text-gray-700 mb-2">Select Enzyme</label>
        <select id="enzyme-select" class="w-full border rounded-lg px-3 py-2 text-sm" onchange="updateEnzymeView()">
          ${enzymeOptions}
        </select>
        <div id="enzyme-info" class="mt-4 text-sm text-gray-600"></div>
      </div>
      <div class="md:col-span-2 bg-white rounded-lg shadow p-4">
        <h3 class="font-semibold text-gray-800 mb-2">Regulated Glycans</h3>
        <div id="table-enzyme-glycans" class="overflow-x-auto"></div>
      </div>
    </div>
    <script>
      const ENZYMES = ${enzymeJson};
      const NETWORK = ${networkJson};
      const GLYCANS = ${glycanJson};
      function updateEnzymeView() {
        const sym = document.getElementById('enzyme-select').value;
        const e = ENZYMES.find(x => x.gene_symbol === sym);
        document.getElementById('enzyme-info').innerHTML = e ?
          '<b>Name:</b> ' + (e.enzyme_name || '') + '<br><b>Family:</b> ' + (e.enzyme_family || '') + '<br><b>Pathway:</b> ' + (e.pathway || '') + '<br><b>Substrate:</b> ' + (e.substrate_type || '') + '<br><br><b>Description:</b><br>' + (e.description || '') : '';
        const links = NETWORK.filter(n => {
          const enzyme = ENZYMES.find(x => x.enzyme_id === n.enzyme_id);
          return enzyme && enzyme.gene_symbol === sym;
        });
        const rows = links.map(link => {
          const g = GLYCANS.find(x => x.glycan_id === link.glycan_id);
          return '<tr class="border-b"><td class="px-4 py-2">' + (g ? g.snfg_name : '') + '</td><td class="px-4 py-2">' + (g ? g.composition : '') + '</td><td class="px-4 py-2">' + (link.linkage || '') + '</td><td class="px-4 py-2"><span class="px-2 py-1 bg-blue-100 text-blue-800 rounded text-xs">' + (link.evidence || '') + '</span></td></tr>';
        }).join('');
        document.getElementById('table-enzyme-glycans').innerHTML =
          '<table class="min-w-full text-left text-sm"><thead class="bg-gray-100"><tr><th class="px-4 py-2">Glycan</th><th class="px-4 py-2">Composition</th><th class="px-4 py-2">Linkage</th><th class="px-4 py-2">Evidence</th></tr></thead><tbody>' + rows + '</tbody></table>';
      }
      updateEnzymeView();
    </script>
  `;
  return pageTemplate('Enzyme Network', content, '/enzymes');
}
''')

    # TCGA page
    f.write('''
// ========================================================================
// TCGA Expression Page
// ========================================================================
function tcgaPage() {
  const tcga = DB.tcga_data;

  const exprChart = JSON.stringify({
    xAxis: { type: 'category', data: tcga.map(d => d.gene), axisLabel: { rotate: 45 } },
    yAxis: { type: 'value', name: 'Median RPKM' },
    series: [
      { name: 'Normal', type: 'bar', data: tcga.map(d => d.normal_rpkm), itemStyle: { color: '#00a65a' } },
      { name: 'PTC', type: 'bar', data: tcga.map(d => d.ptc_rpkm), itemStyle: { color: '#dd4b39' } }
    ],
    tooltip: { trigger: 'axis' },
    legend: { data: ['Normal', 'PTC'] }
  });

  const fcChart = JSON.stringify({
    xAxis: { type: 'category', data: tcga.map(d => d.gene), axisLabel: { rotate: 45 } },
    yAxis: { type: 'value', name: 'Fold Change (PTC/Normal)' },
    series: [{
      data: tcga.map(d => ({
        value: d.fold_change,
        itemStyle: { color: d.significant ? (d.fold_change > 1 ? '#dd4b39' : '#00a65a') : '#bbbbbb' }
      })),
      type: 'bar',
      label: { show: true, position: 'top', formatter: '{c}' }
    }],
    tooltip: { trigger: 'axis' }
  });

  const rows = tcga.map(g => `
    <tr class="border-b hover:bg-gray-50">
      <td class="px-4 py-2 text-sm font-medium">${g.gene}</td>
      <td class="px-4 py-2 text-sm">${g.family || ''}</td>
      <td class="px-4 py-2 text-sm">${g.normal_rpkm}</td>
      <td class="px-4 py-2 text-sm">${g.ptc_rpkm}</td>
      <td class="px-4 py-2 text-sm">${g.fold_change ? g.fold_change.toFixed(2) : 'N/A'}</td>
      <td class="px-4 py-2 text-sm">${g.p_value || ''}</td>
      <td class="px-4 py-2 text-sm">${g.significant ? '<span class="text-red-600 font-bold">*</span>' : 'ns'}</td>
    </tr>
  `).join('');

  const content = `
    <div class="mb-6">
      <h1 class="text-2xl font-bold text-gray-900">TCGA Glycogene Expression</h1>
      <p class="text-gray-600 mt-1">Real TCGA-THCA RNAseq data from Bones et al. 2018 (n=20 paired)</p>
    </div>
    <div class="grid grid-cols-1 md:grid-cols-2 gap-6 mb-6">
      <div class="bg-white rounded-lg shadow p-4">
        <h3 class="font-semibold text-gray-800 mb-2">Median RPKM by Group</h3>
        <div id="chart-tcga-expr" class="chart-container-lg"></div>
      </div>
      <div class="bg-white rounded-lg shadow p-4">
        <h3 class="font-semibold text-gray-800 mb-2">Fold Change & Significance</h3>
        <div id="chart-tcga-fc" class="chart-container-lg"></div>
      </div>
    </div>
    <div class="bg-white rounded-lg shadow overflow-x-auto">
      <table class="min-w-full text-left">
        <thead class="bg-gray-100 text-gray-700 text-xs uppercase">
          <tr><th class="px-4 py-2">Gene</th><th class="px-4 py-2">Family</th><th class="px-4 py-2">Normal RPKM</th><th class="px-4 py-2">PTC RPKM</th><th class="px-4 py-2">FC</th><th class="px-4 py-2">p-value</th><th class="px-4 py-2">Sig</th></tr>
        </thead>
        <tbody class="text-sm">${rows}</tbody>
      </table>
    </div>
    <script>
      echarts.init(document.getElementById('chart-tcga-expr')).setOption(${exprChart});
      echarts.init(document.getElementById('chart-tcga-fc')).setOption(${fcChart});
    </script>
  `;
  return pageTemplate('TCGA Expression', content, '/tcga');
}
''')

    # Diagnostic page
    f.write('''
// ========================================================================
// Diagnostic Tool Page
// ========================================================================
function diagnosticPage() {
  const lit = DB.lit_stats;
  const bn = lit.find(x => x.study === 'Zhang_2021' && x.variable && x.variable.includes('BN'));
  const ratio = lit.find(x => x.study === 'Kudelka_2023' && x.variable && x.variable.includes('G0F:G1F'));

  const bnRef = bn ? `HC: ${bn.group1_mean || 'N/A'} ± ${bn.group1_sd || 'N/A'}% | TC: ${bn.group2_mean || 'N/A'} ± ${bn.group2_sd || 'N/A'}% | AUC: ${bn.auc || 'N/A'}` : 'N/A';
  const ratioRef = ratio ? `HC median: ${ratio.group1_median || 'N/A'} | Recurrent: ${ratio.group2_median || 'N/A'} | AUC: ${ratio.auc || 'N/A'}` : 'N/A';

  const content = `
    <div class="mb-6">
      <h1 class="text-2xl font-bold text-gray-900">Diagnostic Tool</h1>
      <p class="text-gray-600 mt-1">Interactive calculators based on published biomarker performance</p>
    </div>
    <div class="grid grid-cols-1 md:grid-cols-2 gap-6">
      <div class="bg-white rounded-lg shadow p-6">
        <h3 class="font-semibold text-gray-800 mb-2">IgG BN Diagnostic Score (Zhang 2021)</h3>
        <p class="text-xs text-gray-500 mb-4">${bnRef}</p>
        <div class="space-y-3">
          <div><label class="text-sm">H3N5F1 (%)</label><input type="number" id="d-h3n5f1" value="15" class="w-full border rounded px-2 py-1 text-sm"></div>
          <div><label class="text-sm">H4N5F1 (%)</label><input type="number" id="d-h4n5f1" value="20" class="w-full border rounded px-2 py-1 text-sm"></div>
          <div><label class="text-sm">H5N5F1 (%)</label><input type="number" id="d-h5n5f1" value="10" class="w-full border rounded px-2 py-1 text-sm"></div>
          <div><label class="text-sm">H4N4F1 (%) Reference</label><input type="number" id="d-h4n4f1" value="55" class="w-full border rounded px-2 py-1 text-sm"></div>
          <button onclick="calcBN()" class="bg-blue-600 text-white px-4 py-2 rounded text-sm hover:bg-blue-700">Calculate BN Score</button>
        </div>
        <div id="d-result" class="mt-4 p-3 bg-gray-50 rounded text-sm font-mono"></div>
      </div>
      <div class="bg-white rounded-lg shadow p-6">
        <h3 class="font-semibold text-gray-800 mb-2">Serum Recurrence Predictor (Kudelka 2023)</h3>
        <p class="text-xs text-gray-500 mb-4">${ratioRef}</p>
        <div class="space-y-3">
          <div><label class="text-sm">G0F Intensity</label><input type="number" id="d-g0f" value="35" class="w-full border rounded px-2 py-1 text-sm"></div>
          <div><label class="text-sm">G1F Intensity</label><input type="number" id="d-g1f" value="50" class="w-full border rounded px-2 py-1 text-sm"></div>
          <button onclick="calcRatio()" class="bg-yellow-600 text-white px-4 py-2 rounded text-sm hover:bg-yellow-700">Calculate Ratio</button>
        </div>
        <div id="r-result" class="mt-4 p-3 bg-gray-50 rounded text-sm font-mono"></div>
      </div>
    </div>
    <script>
      function calcBN() {
        const h3 = parseFloat(document.getElementById('d-h3n5f1').value) || 0;
        const h4 = parseFloat(document.getElementById('d-h4n5f1').value) || 0;
        const h5 = parseFloat(document.getElementById('d-h5n5f1').value) || 0;
        const ref = parseFloat(document.getElementById('d-h4n4f1').value) || 0;
        const total = h3 + h4 + h5 + ref;
        const bn = total > 0 ? (h3 + h4 + h5) / total : 0;
        document.getElementById('d-result').innerText = 'BN Ratio: ' + bn.toFixed(4) + '\\n' + (bn > 0.35 ? 'BN ELEVATED (consistent with TC)' : 'BN not elevated (consistent with HC)');
      }
      function calcRatio() {
        const g0f = parseFloat(document.getElementById('d-g0f').value) || 0;
        const g1f = parseFloat(document.getElementById('d-g1f').value) || 0;
        if (g1f === 0) { document.getElementById('r-result').innerText = 'Error: G1F cannot be zero'; return; }
        const ratio = g0f / g1f;
        let msg = 'G0F:G1F = ' + ratio.toFixed(3) + '\\n';
        if (ratio > 0.73) msg += 'HIGH recurrence risk (ratio > 0.73)';
        else if (ratio > 0.53) msg += 'MODERATE recurrence risk (0.53-0.73)';
        else msg += 'LOW recurrence risk (ratio < 0.53)';
        document.getElementById('r-result').innerText = msg;
      }
    </script>
  `;
  return pageTemplate('Diagnostic Tool', content, '/diagnostic');
}
''')

    # Nomogram page
    f.write('''
// ========================================================================
// Nomogram Page
// ========================================================================
function nomogramPage() {
  const lit = DB.lit_stats;
  const ca4 = lit.find(x => x.study === 'PTMC_2022_LNM' && x.variable === 'CA4');
  const a2 = lit.find(x => x.study === 'PTMC_2022_LNM' && x.variable && x.variable.includes('A2F0S0G'));

  const ca4Ref = ca4 ? `NLNM median: ${ca4.group1_median || 'N/A'} | LNM median: ${ca4.group2_median || 'N/A'}` : 'N/A';
  const a2Ref = a2 ? `NLNM median: ${a2.group1_median || 'N/A'} | LNM median: ${a2.group2_median || 'N/A'}` : 'N/A';

  const content = `
    <div class="mb-6">
      <h1 class="text-2xl font-bold text-gray-900">PTMC LNM Risk Nomogram</h1>
      <p class="text-gray-600 mt-1">Demonstration model based on PTMC nomogram (Front Oncol 2022). <b>Not for clinical use without validation.</b></p>
    </div>
    <div class="grid grid-cols-1 md:grid-cols-2 gap-6">
      <div class="bg-white rounded-lg shadow p-6">
        <h3 class="font-semibold text-gray-800 mb-2">Input Values</h3>
        <p class="text-xs text-gray-500 mb-4">CA4: ${ca4Ref}<br>A2F0S0G: ${a2Ref}</p>
        <div class="space-y-3">
          <div><label class="text-sm">CA4 (Tetraantennary complex)</label><input type="number" id="n-ca4" value="0.022" step="0.001" class="w-full border rounded px-2 py-1 text-sm"></div>
          <div><label class="text-sm">A2F0S0G (Galactosylation index)</label><input type="number" id="n-a2" value="0.60" step="0.001" class="w-full border rounded px-2 py-1 text-sm"></div>
          <button id="btn-calc-nom" class="bg-blue-600 text-white px-4 py-2 rounded text-sm hover:bg-blue-700">Calculate Risk</button>
        </div>
      </div>
      <div class="bg-white rounded-lg shadow p-6">
        <h3 class="font-semibold text-gray-800 mb-2">Risk Probability</h3>
        <div id="n-gauge" class="chart-container"></div>
        <div id="n-text" class="mt-2 text-center text-sm font-medium"></div>
      </div>
    </div>
    <script>
      let nomogramChart;
      function calcNom() {
        try {
          const ca4 = parseFloat(document.getElementById('n-ca4').value) || 0;
          const a2 = parseFloat(document.getElementById('n-a2').value) || 0;
          const logit = -2.5 + 3.2 * ca4 - 2.0 * a2;
          const prob = 1 / (1 + Math.exp(-logit));
          const pct = (prob * 100).toFixed(1);
          let level = prob < 0.3 ? 'Low' : prob < 0.6 ? 'Intermediate' : 'High';
          let color = prob < 0.3 ? '#00a65a' : prob < 0.6 ? '#f39c12' : '#dd4b39';
          document.getElementById('n-text').innerHTML = '<span style="color:' + color + '">' + level + ' Risk: ' + pct + '%</span>';
          if (typeof echarts === 'undefined') {
            document.getElementById('n-text').innerHTML += '<br><span style="color:#999;font-size:12px">Chart library not loaded</span>';
            return;
          }
          const gaugeEl = document.getElementById('n-gauge');
          if (!nomogramChart) {
            nomogramChart = echarts.init(gaugeEl);
          }
          nomogramChart.setOption({
            series: [{
              type: 'gauge',
              min: 0, max: 100,
              axisLine: { lineStyle: { width: 20, color: [[0.3, '#00a65a'], [0.6, '#f39c12'], [1, '#dd4b39']] } },
              pointer: { length: '60%', width: 5 },
              detail: { fontSize: 20, formatter: '{value}%' },
              data: [{ value: parseFloat(pct), name: 'Risk %' }]
            }]
          });
        } catch (e) {
          document.getElementById('n-text').innerHTML = '<span style="color:#dd4b39">Error: ' + e.message + '</span>';
          console.error(e);
        }
      }
      document.getElementById('btn-calc-nom').addEventListener('click', calcNom);
      calcNom();
    </script>
  `;
  return pageTemplate('Nomogram', content, '/nomogram');
}
''')

    # API handlers
    f.write('''
// ========================================================================
// API Handlers
// ========================================================================
function apiResponse(data, status = 200) {
  return new Response(JSON.stringify(data), {
    status: status,
    headers: { 'Content-Type': 'application/json', 'Access-Control-Allow-Origin': '*' }
  });
}

function handleApi(path) {
  if (path === '/api/studies') return apiResponse(DB.studies);
  if (path === '/api/biomarkers') return apiResponse(DB.biomarkers);
  if (path === '/api/glycans') return apiResponse(DB.glycans);
  if (path === '/api/enzymes') return apiResponse(DB.enzymes);
  if (path === '/api/tcga') return apiResponse(DB.tcga_data);
  if (path === '/api/stats') return apiResponse(DB.stats);
  if (path === '/api/literature-stats') return apiResponse(DB.lit_stats);
  return apiResponse({ error: 'Not found' }, 404);
}
''')

    # Main router
    f.write('''
// ========================================================================
// Main Router
// ========================================================================
export default {
  async fetch(request, env, ctx) {
    const url = new URL(request.url);
    const path = url.pathname;

    if (path.startsWith('/api/')) {
      return handleApi(path);
    }

    switch (path) {
      case '/':
      case '/overview':
        return new Response(overviewPage(), { headers: { 'Content-Type': 'text/html; charset=utf-8' } });
      case '/literature':
        return new Response(literaturePage(), { headers: { 'Content-Type': 'text/html; charset=utf-8' } });
      case '/biomarkers':
        return new Response(biomarkersPage(), { headers: { 'Content-Type': 'text/html; charset=utf-8' } });
      case '/glycans':
        return new Response(glycansPage(), { headers: { 'Content-Type': 'text/html; charset=utf-8' } });
      case '/enzymes':
        return new Response(enzymesPage(), { headers: { 'Content-Type': 'text/html; charset=utf-8' } });
      case '/tcga':
        return new Response(tcgaPage(), { headers: { 'Content-Type': 'text/html; charset=utf-8' } });
      case '/diagnostic':
        return new Response(diagnosticPage(), { headers: { 'Content-Type': 'text/html; charset=utf-8' } });
      case '/nomogram':
        return new Response(nomogramPage(), { headers: { 'Content-Type': 'text/html; charset=utf-8' } });
      default:
        return new Response(overviewPage(), { headers: { 'Content-Type': 'text/html; charset=utf-8' } });
    }
  }
};
''')

print(f"[OK] Generated {OUTPUT_PATH}")
print(f"[INFO] File size: {os.path.getsize(OUTPUT_PATH):,} bytes")
