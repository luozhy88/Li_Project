#!/usr/bin/env python3
"""Regenerate info2-style figures from the project's own dataset and build report v06.

Addresses four revisions:
 1. Top-author table: credit each paper to its first author only (no duplicate articles).
 2. Publication-type breakdown (exclusive categories) + chart, text inserted into section 3.1.
 3. World choropleth map of publications by country + glowing country collaboration network.
 4. Minimalist top-10 countries horizontal bar chart (teal-blue, no grid).
"""
from __future__ import annotations

import csv
import importlib.util
import json
import math
from collections import Counter
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import networkx as nx
import numpy as np
import yaml
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.patches import PathPatch
from matplotlib.path import Path as MplPath
from PIL import Image, ImageDraw, ImageFont
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.shared import Inches, Pt

ROOT = Path('/Users/apple/Desktop/github.local/Li_Project/jiliang')
DATA = ROOT / 'output/01_merged_records.csv'
GEOJSON = ROOT / 'input/mapdata/ne_110m_countries.geojson'
SOURCE = ROOT / 'output/report/05_False_Localizing_Sign_文献计量分析报告_info2图表增强版.docx'
DEST = ROOT / 'output/report/06_False_Localizing_Sign_文献计量分析报告_图表修订版.docx'
ASSET_DIR = ROOT / 'output/figures/geo_style'
PREV_DIR = ROOT / 'output/figures/info2_style'  # reuse author network & journal chart

# Natural Earth NAME differs for a few countries.
NE_NAME = {
    'United States': 'United States of America',
    'Korea, Republic of': 'South Korea',
    'Taiwan, Province of China': 'Taiwan',
}
SHORT = {
    'United States': 'USA', 'United Kingdom': 'UK', 'Korea, Republic of': 'S. KOREA',
    'Taiwan, Province of China': 'TAIWAN', 'Saudi Arabia': 'SAUDI ARABIA',
}


def load_country_stats():
    """Reuse the pipeline's own pycountry-based extraction for consistency."""
    spec = importlib.util.spec_from_file_location('desc', ROOT / '2_descriptive_analysis.py')
    desc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(desc)
    config = yaml.safe_load(open(ROOT / 'config.yaml', encoding='utf-8'))
    patterns = desc.build_country_patterns(config)
    counts, edges = Counter(), Counter()
    for row in read_records():
        cs = sorted(set(desc.extract_countries(row.get('affiliations') or '', patterns)))
        # Merge the two Korea variants produced by alias configuration.
        cs = sorted({'Korea, Republic of' if c == 'South Korea' else c for c in cs})
        counts.update(cs)
        edges.update(combinations(cs, 2))
    return counts, edges


def load_country_year_stats():
    """Per-country publication years (each record counted once per country)."""
    spec = importlib.util.spec_from_file_location('desc', ROOT / '2_descriptive_analysis.py')
    desc = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(desc)
    config = yaml.safe_load(open(ROOT / 'config.yaml', encoding='utf-8'))
    patterns = desc.build_country_patterns(config)
    years = {}
    for row in read_records():
        cs = sorted(set(desc.extract_countries(row.get('affiliations') or '', patterns)))
        cs = sorted({'Korea, Republic of' if c == 'South Korea' else c for c in cs})
        try:
            year = int(float(row.get('year') or 0))
        except (TypeError, ValueError):
            year = 0
        if year <= 0:
            continue
        for c in cs:
            years.setdefault(c, []).append(year)
    return years


def read_records():
    with DATA.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


# ---------------------------------------------------------------- publication types
def classify_types(records):
    """Exclusive single-label classification by priority."""
    groups = [
        ('Case Reports', {'Case Reports'}),
        ('Review (incl. Systematic Review)', {'Review', 'Systematic Review'}),
        ('Clinical Trial', {'Clinical Trial', 'Randomized Controlled Trial'}),
        ('Letter / Editorial / Comment', {'Letter', 'Editorial', 'Comment'}),
        ('Historical / Biography', {'Historical Article', 'Biography', 'Portrait'}),
    ]
    fallback = {'Journal Article', 'J', 'English Abstract', 'Comparative Study',
                'Research Support, Non-U.S. Gov\'t', 'Review Journal Article'}
    counts = Counter()
    for r in records:
        types = {t.strip() for t in (r.get('publication_types') or '').split(';') if t.strip()}
        label = None
        for name, keys in groups:
            if types & keys:
                label = name
                break
        if label is None:
            label = 'Journal Article (other)' if (types & fallback or not types) else 'Other'
        counts[label] += 1
    order = [g[0] for g in groups] + ['Journal Article (other)', 'Other']
    return [(k, counts[k]) for k in order if counts[k]]


# ---------------------------------------------------------------- fig E1: world choropleth
def draw_world_map(path, counts):
    geo = json.load(open(GEOJSON, encoding='utf-8'))
    value_by_ne = {}
    for name, n in counts.items():
        value_by_ne[NE_NAME.get(name, name)] = n
    unmatched = [k for k in value_by_ne
                 if not any(f['properties'].get('NAME') == k for f in geo['features'])]
    if unmatched:
        print('  [warn] no map feature for:', unmatched)
    vmax = max(counts.values())
    cmap = plt.get_cmap('Blues')
    norm = Normalize(vmin=0, vmax=vmax)

    fig, ax = plt.subplots(figsize=(14.5, 7.6), dpi=200)
    fig.patch.set_facecolor('white')
    centroids = {}
    for feat in geo['features']:
        name = feat['properties'].get('NAME')
        geom = feat.get('geometry') or {}
        polys = geom.get('coordinates') or []
        if geom.get('type') == 'Polygon':
            polys = [polys]
        value = value_by_ne.get(name)
        color = '#d9d9d9' if value is None else cmap(norm(value))
        best_area, best_ring = -1.0, None
        for poly in polys:
            ring = poly[0]
            pts = [(x, y) for x, y in ring]
            verts = pts + [pts[0]]
            codes = [MplPath.MOVETO] + [MplPath.LINETO] * (len(pts) - 1) + [MplPath.CLOSEPOLY]
            ax.add_patch(PathPatch(MplPath(verts, codes), facecolor=color,
                                   edgecolor='white', linewidth=.4, zorder=2))
            area = abs(sum(pts[i][0] * pts[(i + 1) % len(pts)][1]
                           - pts[(i + 1) % len(pts)][0] * pts[i][1]
                           for i in range(len(pts))) / 2)
            if area > best_area:
                best_area, best_ring = area, pts
        if value is not None and best_ring:
            a = sum(best_ring[i][0] * best_ring[(i + 1) % len(best_ring)][1]
                    - best_ring[(i + 1) % len(best_ring)][0] * best_ring[i][1]
                    for i in range(len(best_ring))) / 2
            cx = sum((best_ring[i][0] + best_ring[(i + 1) % len(best_ring)][0])
                     * (best_ring[i][0] * best_ring[(i + 1) % len(best_ring)][1]
                        - best_ring[(i + 1) % len(best_ring)][0] * best_ring[i][1])
                     for i in range(len(best_ring))) / (6 * a)
            cy = sum((best_ring[i][1] + best_ring[(i + 1) % len(best_ring)][1])
                     * (best_ring[i][0] * best_ring[(i + 1) % len(best_ring)][1]
                        - best_ring[(i + 1) % len(best_ring)][0] * best_ring[i][1])
                     for i in range(len(best_ring))) / (6 * a)
            centroids[name] = (cx, cy)

    offsets = {  # fine-tune crowded labels (lon, lat shifts)
        'United Kingdom': (-2, 3.5), 'Netherlands': (0, 2.5), 'Switzerland': (-2.5, 2.5),
        'Slovakia': (2.5, 1.5), 'Denmark': (3, 1.5), 'Israel': (1.5, 2.5),
        'South Korea': (3.5, -3), 'Taiwan': (3, 0), 'Malaysia': (2, -4),
        'Hungary': (0, -2.5), 'Georgia': (4, 1), 'Italy': (1.5, -1), 'Nepal': (0.5, 1.5),
    }
    for name, value in value_by_ne.items():
        if name not in centroids:
            continue
        x, y = centroids[name]
        dx, dy = offsets.get(name, (0, 0))
        shade = norm(value)
        ax.text(x + dx, y + dy, str(value), ha='center', va='center', zorder=3,
                fontsize=7.5 if value < 10 else 8.5, fontweight='bold',
                color='white' if shade > .45 else '#08306b')
    ax.set_xlim(-170, 195); ax.set_ylim(-58, 84)
    ax.axis('off')

    cax = fig.add_axes([0.055, 0.22, 0.014, 0.5])
    sm = ScalarMappable(norm=Normalize(vmin=1, vmax=vmax), cmap=cmap)
    cb = fig.colorbar(sm, cax=cax)
    cb.set_ticks([1, round(vmax / 2), vmax])
    cb.ax.tick_params(labelsize=9)
    cb.set_label('Freq', fontsize=11)
    cb.outline.set_visible(False)
    fig.subplots_adjust(left=.01, right=.99, top=.99, bottom=.01)
    fig.savefig(path, facecolor='white', bbox_inches='tight')
    plt.close(fig)


# ---------------------------------------------------------------- fig E2: collaboration network
def draw_glow_network(path, counts, edges, country_years=None):
    """Country collaboration network with CiteSpace-style annual-ring nodes.

    - Node AREA strictly proportional to total publication output.
    - Each node is drawn as concentric rings: one ring per decade with at
      least one publication, ring color encodes the decade (shared color
      scale), ring area is proportional to the publication count in that
      decade.  Earliest decade = innermost ring, latest = outermost.
    - Connected components are hand-placed to avoid label overlap and edge
      crossings; isolated nodes are aligned in a bottom strip.
    - A legend (bottom right) explains node size, ring colors and links.
    """
    g = nx.Graph()
    for (a, b), w in edges.items():
        g.add_edge(a, b, weight=w)
    connected = set(g.nodes)
    isolated = sorted((set(counts) - connected), key=lambda c: (-counts[c], c))

    # --- hand-placed core layout: one quadrant per connected component,
    #     edges fan out from each hub so no two links cross
    pos = {
        # USA component (upper centre-left)
        'United States': (-0.30, 1.08),
        'Mexico': (-1.10, 1.42),
        'Italy': (0.88, 1.32),
        'Slovakia': (-0.22, 0.52),
        # UK component (upper right)
        'United Kingdom': (1.28, 0.98),
        'Australia': (1.74, 1.40),
        'Denmark': (1.76, 0.56),
        # Germany--Switzerland (lower left)
        'Germany': (-1.12, 0.28),
        'Switzerland': (-1.58, 0.02),
        # S. Korea--Malaysia (lower centre-right)
        'Korea, Republic of': (0.52, 0.34),
        'Malaysia': (1.04, 0.06),
    }
    # --- isolated nodes: bottom strip, two rows sorted by output
    row1, row2 = isolated[:8], isolated[8:]
    for i, name in enumerate(row1):
        pos[name] = (-1.62 + 2.50 * i / max(len(row1) - 1, 1), -0.82)
    for i, name in enumerate(row2):
        pos[name] = (-1.50 + 2.38 * i / max(len(row2) - 1, 1), -1.34)

    fig, ax = plt.subplots(figsize=(13.5, 9), dpi=200)
    fig.patch.set_facecolor('white')
    vmax = max(counts.values())

    def size(c):  # scatter area in points^2, strictly linear in c
        return 100 + 210 * c

    cmap = plt.get_cmap('turbo')

    def ring_periods(name):
        """{5-year slice: count} for a country, from its publication years."""
        ys = (country_years or {}).get(name)
        if not ys:  # fallback: single ring of neutral colour
            return {2025: counts[name]}
        pc = Counter(y // 5 * 5 for y in ys)
        return dict(sorted(pc.items()))

    # color scale normalised to the slices actually present in the data
    node_periods = {name: ring_periods(name) for name in counts}
    all_slices = [s for periods in node_periods.values() for s in periods]
    smin, smax = min(all_slices), max(all_slices)
    norm = Normalize(smin, smax)

    # links
    for (a, b), w in edges.items():
        x1, y1 = pos[a]; x2, y2 = pos[b]
        ax.plot([x1, x2], [y1, y2], color='#d73027', alpha=.5,
                linewidth=.9 + .8 * w, zorder=1)

    # nodes: concentric annual rings (outermost = latest decade)
    for name in counts:
        x, y = pos[name]
        s_total = size(counts[name])
        periods = node_periods[name]
        total = sum(periods.values())
        cum = {}
        run = 0
        for dec in sorted(periods):
            run += periods[dec]
            cum[dec] = run / total
        for dec in sorted(periods, reverse=True):  # largest disc first
            ax.scatter([x], [y], s=s_total * cum[dec], color=cmap(norm(dec)),
                       linewidths=.5, edgecolors='white', zorder=3)
        label = SHORT.get(name, name).upper()
        ax.text(x, y, label, ha='center', va='center', zorder=6,
                fontsize=5.5 + 4.5 * math.sqrt(counts[name]) / math.sqrt(vmax),
                fontweight='bold', color='#2b2b2b',
                path_effects=[pe.withStroke(linewidth=1.6, foreground='white')])

    # divider between core network and isolated-node strip
    ax.plot([-1.95, 2.25], [-0.48, -0.48], color='#999999', lw=.8,
            ls=(0, (4, 3)), zorder=1)
    ax.text(-1.93, -0.56, 'No co-authorship links (isolated nodes)',
            fontsize=8.5, color='#666666', va='center', style='italic')

    # --- legend (bottom right)
    lx, ly = 1.10, -1.86  # box lower-left corner
    ax.add_patch(plt.Rectangle((lx, ly), 1.55, 1.60, facecolor='white',
                               edgecolor='#999999', lw=.9, zorder=7))
    ax.text(lx + .775, ly + 1.48, 'Legend', fontsize=10, fontweight='bold',
            ha='center', zorder=8)
    ax.text(lx + .10, ly + 1.28, 'Node area ∝ total publication output:',
            fontsize=8, va='center', zorder=8)
    for cx, c in [(lx + .42, 1), (lx + .80, 8), (lx + 1.22, 24)]:
        ax.scatter([cx], [ly + 1.02], s=size(c), color='#8c8c8c', alpha=.9,
                   linewidths=.4, edgecolors='white', zorder=8)
        ax.text(cx, ly + .70, str(c), fontsize=8, ha='center', zorder=8)
    ax.text(lx + .10, ly + .54, 'Ring color = publication period:',
            fontsize=8, va='center', zorder=8)
    gradient = np.linspace(0, 1, 256).reshape(1, -1)
    ax.imshow(gradient, cmap=cmap, aspect='auto', zorder=8,
              extent=(lx + .10, lx + 1.45, ly + .36, ly + .44))
    ax.add_patch(plt.Rectangle((lx + .10, ly + .36), 1.35, .08, fill=False,
                               edgecolor='#999999', lw=.6, zorder=9))
    ax.text(lx + .10, ly + .24, str(smin), fontsize=7.5, ha='left',
            va='center', zorder=8)
    ax.text(lx + 1.45, ly + .24, str(smax), fontsize=7.5, ha='right',
            va='center', zorder=8)
    ax.plot([lx + .10, lx + .32], [ly + .10, ly + .10], color='#d73027',
            alpha=.6, lw=1.7, zorder=8)
    ax.text(lx + .38, ly + .10, 'Co-authorship link (shared record)',
            fontsize=8, va='center', zorder=8)

    ax.text(.01, .985, 'Network visualization of international collaboration among countries',
            transform=ax.transAxes, fontsize=14, fontweight='bold', va='top')
    ax.text(.01, .945, 'Node area ∝ publication output; ring color = period; '
            'link = co-authorship.',
            transform=ax.transAxes, fontsize=9.5, color='#444444', va='top')
    ax.axis('off')
    ax.set_xlim(-1.98, 2.78)
    ax.set_ylim(-1.96, 1.72)
    fig.savefig(path, facecolor='white', bbox_inches='tight')
    plt.close(fig)


# ---------------------------------------------------------------- fig E3: top-10 countries bar
def draw_top10_countries(path, counts):
    top = counts.most_common(10)
    fig, ax = plt.subplots(figsize=(10.5, 6.2), dpi=200)
    fig.patch.set_facecolor('white')
    names = [SHORT.get(c, c) for c, _ in top][::-1]
    vals = [n for _, n in top][::-1]
    ax.barh(names, vals, color='#0f7e8c', height=.62)
    ax.set_title('Top 10 countries/regions publishing on false localizing signs',
                 fontsize=13, fontweight='bold', pad=12)
    ax.set_xlabel('Number of publications', fontsize=11)
    ax.set_ylabel('Country / Region', fontsize=11)
    ax.set_xlim(0, 25)
    ax.set_xticks(range(0, 26, 5))
    ax.grid(False)
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    for i, v in enumerate(vals):
        ax.text(v + .35, i, str(v), va='center', fontsize=10, color='#0f7e8c',
                fontweight='bold')
    fig.tight_layout()
    fig.savefig(path, facecolor='white', bbox_inches='tight')
    plt.close(fig)


# ---------------------------------------------------------------- fig E4: publication types
def draw_type_chart(path, type_counts, total):
    fig, ax = plt.subplots(figsize=(10.5, 5.6), dpi=200)
    fig.patch.set_facecolor('white')
    labels = [k for k, _ in type_counts][::-1]
    vals = [n for _, n in type_counts][::-1]
    colors = ['#4d9ed7', '#e58d47', '#47a85c', '#a46ccc', '#d2af34', '#e85d5d', '#8c9aa5'][::-1]
    ax.barh(labels, vals, color=colors[:len(vals)], height=.6)
    ax.set_title('Publication types of the included records', fontsize=13,
                 fontweight='bold', pad=12)
    ax.set_xlabel('Number of publications', fontsize=11)
    ax.grid(False)
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    for i, v in enumerate(vals):
        ax.text(v + .8, i, f'{v}  ({v / total * 100:.1f}%)', va='center',
                fontsize=10, color='#333333')
    ax.set_xlim(0, max(vals) * 1.22)
    fig.tight_layout()
    fig.savefig(path, facecolor='white', bbox_inches='tight')
    plt.close(fig)


# ---------------------------------------------------------------- table E1: top first authors
def get_font(size, bold=False):
    for item in ['/System/Library/Fonts/Supplemental/Arial Bold.ttf' if bold else
                 '/System/Library/Fonts/Supplemental/Arial.ttf',
                 '/System/Library/Fonts/Supplemental/Arial.ttf']:
        if Path(item).exists():
            return ImageFont.truetype(item, size)
    return ImageFont.load_default()


def wrap(draw, text, font, width):
    words, out, line = text.split(), [], ''
    for word in words:
        test = (line + ' ' + word).strip()
        if line and draw.textlength(test, font=font) > width:
            out.append(line)
            line = word
        else:
            line = test
    return out + ([line] if line else [])


def draw_first_author_table(path, records, label='Table D1'):
    counts = Counter(r.get('first_author', '').strip() for r in records
                     if r.get('first_author', '').strip())
    top = counts.most_common(10)
    W, H = 1800, 1320
    im = Image.new('RGB', (W, H), '#fff9bd')
    draw = ImageDraw.Draw(im)
    title, head, body, tiny = get_font(34, True), get_font(25, True), get_font(23), get_font(17)
    draw.text((60, 42), f'{label}. Top 10 productive first authors in false localizing sign research',
              font=title, fill='#111')
    draw.line((60, 108, W - 60, 108), fill='#111', width=3)
    for x, lab in zip([60, 420, 690],
                      ['First author', 'Documents', 'Representative first-authored publication']):
        draw.text((x, 128), lab, font=head, fill='#111')
    draw.line((60, 174, W - 60, 174), fill='#111', width=2)
    y = 192
    for author, n in top:
        own = [r for r in records if r.get('first_author', '').strip() == author]
        chosen = sorted(own, key=lambda r: int(r.get('year') or 0), reverse=True)[0]
        article = ' '.join((chosen.get('title') or '').split())
        lines = wrap(draw, article, body, 1030)[:3]
        row_h = max(82, 25 + 30 * len(lines))
        draw.text((60, y + 12), author, font=body, fill='#111')
        draw.text((430, y + 12), f'{n} ({n / len(records) * 100:.2f}%)', font=body, fill='#111')
        for j, ln in enumerate(lines):
            draw.text((690, y + 12 + j * 30), ln, font=body, fill='#111')
        draw.line((60, y + row_h, W - 60, y + row_h), fill='#777', width=1)
        y += row_h
    draw.text((60, H - 60),
              f'Each of the {len(records)} records is credited to its first author only, so every '
              'publication appears exactly once; the final column lists one recent first-authored '
              'paper per author.', font=tiny, fill='#333')
    im.save(path)
    return top


# ---------------------------------------------------------------- docx assembly
def set_font(run):
    run.font.name = 'Songti SC'
    run._element.rPr.rFonts.set(qn('w:eastAsia'), 'Songti SC')


def add_text(doc, text):
    p = doc.add_paragraph()
    set_font(p.add_run(text))
    return p


def add_heading(doc, text, level):
    p = doc.add_heading(level=level)
    set_font(p.add_run(text))
    return p


def add_image(doc, path, caption):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.add_run().add_picture(str(path), width=Inches(6.25))
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run(caption)
    run.italic = True
    run.font.size = Pt(9)
    set_font(run)


def insert_image_after(doc, anchor_par, path, caption):
    """Insert an image paragraph + caption right after an existing paragraph."""
    cap = doc.add_paragraph()
    cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = cap.add_run(caption)
    run.italic = True
    run.font.size = Pt(9)
    set_font(run)
    pic = doc.add_paragraph()
    pic.alignment = WD_ALIGN_PARAGRAPH.CENTER
    pic.add_run().add_picture(str(path), width=Inches(6.0))
    anchor_par._p.addnext(cap._p)
    anchor_par._p.addnext(pic._p)


def main():
    records = read_records()
    ASSET_DIR.mkdir(parents=True, exist_ok=True)
    counts, edges = load_country_stats()
    type_counts = classify_types(records)
    total = len(records)
    print('文献类型（互斥分类）:', type_counts)
    print('国家数:', len(counts), '| Top:', counts.most_common(5))

    assets = {
        'map': ASSET_DIR / 'fig_e1_world_map.png',
        'network': ASSET_DIR / 'fig_e2_country_network_glow.png',
        'top10': ASSET_DIR / 'fig_e3_top10_countries.png',
        'types': ASSET_DIR / 'fig_e4_publication_types.png',
        'table': ASSET_DIR / 'table_e1_top_first_authors.png',
    }
    draw_world_map(assets['map'], counts)
    draw_glow_network(assets['network'], counts, edges,
                      country_years=load_country_year_stats())
    draw_top10_countries(assets['top10'], counts)
    draw_type_chart(assets['types'], type_counts, total)
    top_first = draw_first_author_table(assets['table'], records)
    print('Top first authors:', top_first)

    doc = Document(SOURCE)

    # ---- issue 2: publication-type sentence + chart in section 3.1
    zh = {'Case Reports': '病例报告', 'Review (incl. Systematic Review)': '综述（含系统综述）',
          'Clinical Trial': '临床试验', 'Letter / Editorial / Comment': '来信/述评/评论',
          'Historical / Biography': '历史/传记类', 'Journal Article (other)': '其他期刊论文',
          'Other': '其他类型'}
    parts = '，'.join(f'{zh.get(k, k)}{n} 篇（{n / total * 100:.1f}%）' for k, n in type_counts)
    sentence = f'文献体裁构成方面（按优先级互斥归类）：{parts}。'
    anchor = None
    for p in doc.paragraphs:
        if '具有摘要的记录' in p.text:
            anchor = p
            break
    assert anchor is not None, '未找到 3.1 节目标段落'
    run = anchor.add_run(sentence)
    set_font(run)
    insert_image_after(doc, anchor, assets['types'],
                       'Fig. 3-1. Distribution of publication types among the included records '
                       '(exclusive single-label classification).')

    # ---- remove the old appendix D (keep the final sectPr)
    body = doc.element.body
    appendix_el = None
    for p in doc.paragraphs:
        if p.text.strip().startswith('附录 D'):
            appendix_el = p._p
            break
    assert appendix_el is not None, '未找到附录 D 标题'
    for sib in list(appendix_el.itersiblings()):
        if sib.tag == qn('w:sectPr'):
            continue
        body.remove(sib)
    body.remove(appendix_el)

    # ---- rebuild appendix D with the new figures
    doc.add_page_break()
    add_heading(doc, '附录 D  国家/地区分布、合作网络与核心作者分析（修订版）', 1)
    add_text(doc, '本附录参照示例文章中国家分布地图、国家合作网络、核心作者表和期刊条形图的版式，'
                  '使用本研究纳入的 159 篇假定位体征文献重新绘制。所有统计均基于合并去重后的题录数据；'
                  '国家/地区根据作者地址字段识别，合作连线仅表示同一篇题录中的共同署名，'
                  '不代表引文关系或机构协作强度。')

    add_heading(doc, 'D.1 国家/地区发文分布', 2)
    add_text(doc, f'图 D1 以单色调蓝色渐变的世界地图展示各国/地区发文量，颜色越深表示发文越多'
                  f'（灰色为无数据）；图 D2 以横向条形图展示发文量前 10 位的国家/地区。'
                  f'美国以 {counts["United States"]} 篇居首，其次为印度（{counts["India"]} 篇）、'
                  f'日本（{counts["Japan"]} 篇）和英国（{counts["United Kingdom"]} 篇）。')
    add_image(doc, assets['map'],
              'Fig. D1. Geographic distribution of false localizing sign publications worldwide. '
              'The color intensity reflects the publication frequency of each country, with darker '
              'shades indicating higher research productivity (grey: no records).')
    doc.add_page_break()
    add_image(doc, assets['top10'],
              'Fig. D2. Top 10 countries/regions publishing on false localizing signs.')

    doc.add_page_break()
    add_heading(doc, 'D.2 国家/地区合作网络', 2)
    add_text(doc, '图 D3 展示国家/地区间的合作网络：节点大小与发文量成正比，节点采用黄—橙—红多层配色，'
                  '红色连线表示两国/地区在同一篇文献中共同署名。')
    add_image(doc, assets['network'],
              'Fig. D3. Network visualization of international collaboration among countries in '
              'false localizing sign research. Each node represents a country, with node size '
              'proportional to publication output; links indicate collaborative relationships '
              'between countries.')

    doc.add_page_break()
    add_heading(doc, 'D.3 核心第一作者与作者合作网络', 2)
    add_text(doc, '为避免同一篇文献在多位作者间重复计数，表 D1 按第一作者统计发文量，'
                  '每篇文献仅计入其第一作者一次；图 D4 展示存在共同署名关系的作者合作核心。')
    add_image(doc, assets['table'],
              'Table D1. Top 10 productive first authors in false localizing sign research.')
    doc.add_page_break()
    add_image(doc, PREV_DIR / 'fig_d2_author_collaboration_network.png',
              'Fig. D4. Co-authorship network of authors in false localizing sign research.')

    doc.add_page_break()
    add_heading(doc, 'D.4 核心期刊', 2)
    add_text(doc, '图 D5 展示本数据集发文量前 15 位的期刊。')
    add_image(doc, PREV_DIR / 'fig_d3_top_journals.png',
              'Fig. D5. Top 15 journals publishing on false localizing signs.')

    doc.save(DEST)
    print(DEST)


if __name__ == '__main__':
    main()
