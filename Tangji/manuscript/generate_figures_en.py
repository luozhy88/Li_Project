#!/usr/bin/env python3
"""
generate_figures_en.py
Generate English-labeled figures for SCI manuscript
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import numpy as np
import sqlite3
import os

DB_PATH = "output/thyroid_glyco_db.sqlite"
FIG_DIR = "manuscript/figures"
os.makedirs(FIG_DIR, exist_ok=True)

con = sqlite3.connect(DB_PATH)

# ========================================================================
# Figure 1: Literature trend
# ========================================================================
def fig1():
    cursor = con.cursor()
    cursor.execute("SELECT year, COUNT(*) FROM studies WHERE year IS NOT NULL GROUP BY year ORDER BY year")
    rows = cursor.fetchall()
    years = [r[0] for r in rows]
    counts = [r[1] for r in rows]
    
    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(years, counts, color='#3c8dbc', edgecolor='white', width=0.6)
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.05,
                f'{int(count)}', ha='center', va='bottom', fontsize=10)
    ax.set_xlabel('Publication Year', fontsize=12)
    ax.set_ylabel('Number of Publications', fontsize=12)
    ax.set_title('Figure 1  Annual publications on thyroid cancer glycosylation', fontsize=13, fontweight='bold')
    ax.set_xticks(years)
    ax.set_ylim(0, max(counts) + 1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/fig1_literature_trend_en.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("[OK] Figure 1 generated")

# ========================================================================
# Figure 2: Biomarker AUC comparison
# ========================================================================
def fig2():
    cursor = con.cursor()
    cursor.execute("SELECT name, performance_auc FROM biomarkers WHERE performance_auc IS NOT NULL")
    rows = cursor.fetchall()
    names = [r[0] for r in rows]
    aucs = [r[1] for r in rows]
    
    fig, ax = plt.subplots(figsize=(8, 5))
    colors = ['#dd4b39' if a >= 0.9 else '#f39c12' if a >= 0.8 else '#999999' for a in aucs]
    bars = ax.barh(range(len(names)), aucs, color=colors, edgecolor='white', height=0.5)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel('AUC', fontsize=12)
    ax.set_title('Figure 2  Diagnostic performance of reported glycosylation biomarkers', fontsize=13, fontweight='bold')
    ax.axvline(x=0.8, color='orange', linestyle='--', linewidth=1, label='AUC=0.80')
    ax.axvline(x=0.9, color='green', linestyle='--', linewidth=1, label='AUC=0.90')
    ax.set_xlim(0.7, 1.0)
    for i, (bar, auc) in enumerate(zip(bars, aucs)):
        ax.text(auc + 0.01, i, f'{auc:.3f}', va='center', fontsize=9)
    ax.legend(loc='lower right', fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/fig2_biomarker_auc_en.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("[OK] Figure 2 generated")

# ========================================================================
# Figure 3: TCGA glycogene expression
# ========================================================================
def fig3():
    cursor = con.cursor()
    cursor.execute("SELECT gene_symbol, normal_median_rpkm, ptc_median_rpkm, significant FROM tcga_glycogene_expression ORDER BY fold_change DESC")
    rows = cursor.fetchall()
    genes = [r[0] for r in rows]
    normal = [r[1] for r in rows]
    ptc = [r[2] for r in rows]
    sig = [r[3] for r in rows]
    
    x = np.arange(len(genes))
    width = 0.35
    fig, ax = plt.subplots(figsize=(10, 6))
    bars1 = ax.bar(x - width/2, normal, width, label='Normal (n=20)', color='#00a65a', edgecolor='white')
    bars2 = ax.bar(x + width/2, ptc, width, label='PTC (n=20)', color='#dd4b39', edgecolor='white')
    for i, s in enumerate(sig):
        if s == 1:
            max_h = max(normal[i], ptc[i])
            ax.plot([i-width/2, i+width/2], [max_h+0.3, max_h+0.3], 'k-', linewidth=1)
            ax.text(i, max_h+0.45, '*', ha='center', fontsize=14)
    ax.set_ylabel('Median RPKM', fontsize=12)
    ax.set_title('Figure 3  Differential expression of glycosyltransferases in TCGA-THCA\n(Source: Bones et al. Cancers 2018)', fontsize=12, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(genes, rotation=45, ha='right', fontsize=10)
    ax.legend(loc='upper right', fontsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/fig3_tcga_glycogene_en.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("[OK] Figure 3 generated")

# ========================================================================
# Figure 4: Platform architecture
# ========================================================================
def fig4():
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 10)
    ax.axis('off')
    ax.text(6, 9.5, 'Figure 4  ThyGlycoPortal platform architecture', fontsize=15, fontweight='bold', ha='center')
    
    c_data, c_ana, c_cli, c_user = '#3498db', '#e74c3c', '#2ecc71', '#f39c12'
    
    def draw_box(x, y, w, h, text, color, fs=9):
        box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05", facecolor=color, edgecolor='black', alpha=0.8)
        ax.add_patch(box)
        ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=fs, color='white', fontweight='bold')
    
    draw_box(0.5, 7.5, 2.5, 1.2, 'Literature DB\n(11 studies)', c_data)
    draw_box(3.5, 7.5, 2.5, 1.2, 'Glycan Library\n(24 N-glycans)', c_data)
    draw_box(6.5, 7.5, 2.5, 1.2, 'TCGA Expression\n(20 paired)', c_data)
    draw_box(9.5, 7.5, 2.0, 1.2, 'User Upload\nData', c_data)
    ax.text(6, 7.0, 'Data Resource Layer', fontsize=11, ha='center', fontweight='bold')
    
    draw_box(0.5, 4.8, 2.0, 1.2, 'Glycan\nBrowser', c_ana)
    draw_box(3.0, 4.8, 2.0, 1.2, 'Biomarker\nValidation', c_ana)
    draw_box(5.5, 4.8, 2.0, 1.2, 'Enzyme-Glycan\nNetwork', c_ana)
    draw_box(8.0, 4.8, 2.0, 1.2, 'Multi-omics\nIntegration', c_ana)
    draw_box(10.2, 4.8, 1.5, 1.2, 'Meta\nAnalysis', c_ana)
    ax.text(6, 4.3, 'Core Analysis Layer', fontsize=11, ha='center', fontweight='bold')
    
    draw_box(2.0, 2.5, 3.0, 1.2, 'Nomogram Tool\n(LNM Risk)', c_cli)
    draw_box(6.0, 2.5, 3.0, 1.2, 'Clinical Decision\nSupport', c_cli)
    ax.text(6, 2.0, 'Clinical Decision Layer', fontsize=11, ha='center', fontweight='bold')
    
    draw_box(4.0, 0.5, 4.0, 1.0, 'End-user Interface (R Shiny)', c_user, fs=11)
    
    arrow = dict(arrowstyle='->', color='gray', lw=1.5)
    for x in [1.75, 4.75, 7.75, 10.5]:
        ax.annotate('', xy=(x, 6.0), xytext=(x, 7.5), arrowprops=arrow)
    for x in [1.5, 4.0, 6.5, 9.0, 10.95]:
        ax.annotate('', xy=(x, 3.7), xytext=(x, 4.8), arrowprops=arrow)
    ax.annotate('', xy=(6, 1.5), xytext=(6, 2.5), arrowprops=arrow)
    
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/fig4_platform_architecture_en.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("[OK] Figure 4 generated")

# ========================================================================
# Figure 5: Glycan change patterns
# ========================================================================
def fig5():
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    
    ax = axes[0]
    categories = ['HC\n(n=15)', 'Recurrent DTC\n(n=13)']
    g0f = [22.5, 31.2]
    g1f = [38.5, 29.8]
    g2f = [18.2, 13.5]
    x = np.arange(len(categories))
    width = 0.25
    ax.bar(x - width, g0f, width, label='G0F', color='#e74c3c')
    ax.bar(x, g1f, width, label='G1F', color='#3498db')
    ax.bar(x + width, g2f, width, label='G2F', color='#2ecc71')
    ax.set_ylabel('Relative abundance (%)', fontsize=10)
    ax.set_title('A. Recurrence Prediction\n(Kudelka 2023)', fontsize=11, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=9)
    ax.legend(fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    ax = axes[1]
    categories = ['HC\n(n=25)', 'TC\n(n=25)']
    bn = [18.5, 27.5]
    ax.bar(categories, bn, color=['#00a65a', '#dd4b39'], edgecolor='white', width=0.5)
    ax.set_ylabel('BN relative abundance (%)', fontsize=10)
    ax.set_title('B. IgG BN Diagnosis\n(Zhang 2021, AUC=0.920)', fontsize=11, fontweight='bold')
    for i, v in enumerate(bn):
        ax.text(i, v + 0.5, f'{v}%', ha='center', fontsize=10, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    ax = axes[2]
    genes = ['ST6GAL1', 'MAN1A2', 'MAN2A1', 'B4GALT1', 'MGAT5']
    fc = [1.62, 0.69, 0.73, 0.74, 1.50]
    colors = ['#dd4b39' if f > 1 else '#00a65a' for f in fc]
    bars = ax.barh(range(len(genes)), fc, color=colors, edgecolor='white')
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=10)
    ax.set_xlabel('Fold change (PTC/Normal)', fontsize=10)
    ax.set_title('C. Glycosyltransferase Alterations\n(Bones 2018, TCGA)', fontsize=11, fontweight='bold')
    ax.axvline(x=1, color='black', linestyle='--', linewidth=1)
    for i, (bar, f) in enumerate(zip(bars, fc)):
        ax.text(f + 0.05 if f > 1 else f - 0.05, i, f'{f:.2f}', va='center', ha='left' if f > 1 else 'right', fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.suptitle('Figure 5  Summary of glycosylation alterations in thyroid cancer', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/fig5_glycan_changes_en.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("[OK] Figure 5 generated")

# ========================================================================
# Table 1: Database comparison (as figure)
# ========================================================================
def table1():
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.axis('tight')
    ax.axis('off')
    data = [
        ['Database', 'Type', 'Thyroid-specific', 'Glycosylation', 'Access'],
        ['GlycoPOST', 'MS repository', 'No', 'Raw MS', 'Free'],
        ['UniCarbKB', 'Glycan KB', 'No', 'Structures', 'Free'],
        ['GlycoStore', 'LC-MS DB', 'No', 'Experimental', 'Free'],
        ['TCGA', 'Cancer multi-omics', 'Partial (genomics)', 'Minimal', 'Free'],
        ['THPA', 'Protein atlas', 'Partial (IHC)', 'No', 'Free'],
        ['GlyConnect', 'Glycoproteomics', 'No', 'Glycoproteins', 'Free'],
        ['ThyGlycoPortal (This work)', 'Thyroid glyco-specific', 'Yes', 'Integrated analysis', 'Free']
    ]
    table = ax.table(cellText=data, cellLoc='center', loc='center', colWidths=[0.25, 0.20, 0.20, 0.20, 0.15])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2)
    for i in range(5):
        cell = table[(0, i)]
        cell.set_facecolor('#3c8dbc')
        cell.set_text_props(weight='bold', color='white')
    for i in range(5):
        cell = table[(7, i)]
        cell.set_facecolor('#d4edda')
        cell.set_text_props(weight='bold')
    ax.set_title('Table 1  Comparison of existing databases with ThyGlycoPortal', fontsize=13, fontweight='bold', pad=20)
    plt.savefig(f'{FIG_DIR}/table1_database_comparison_en.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("[OK] Table 1 generated")

# ========================================================================
# Table 2: Biomarker summary (as figure)
# ========================================================================
def table2():
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.axis('tight')
    ax.axis('off')
    data = [
        ['Biomarker', 'Sample', 'Comparison', 'AUC', 'p-value', 'Key Finding'],
        ['BN (bisecting neutral)', 'Plasma IgG', 'TC vs HC/BTN', '0.920', '<0.0001', 'Elevated in TC'],
        ['Glyco-panel', 'Plasma IgG', 'TC vs HC/BTN', '0.917', '<0.0001', '4-glycan classifier'],
        ['G0F:G1F ratio', 'Serum', 'Recurrent vs HC', '0.820', '0.004', 'Hypogalactosylation predicts recurrence'],
        ['CA4 + A2F0S0G', 'Serum', 'LNM vs NLNM', '0.702/0.658', '0.001/0.011', 'Branching/galactosylation predict LNM']
    ]
    table = ax.table(cellText=data, cellLoc='center', loc='center', colWidths=[0.20, 0.12, 0.18, 0.10, 0.12, 0.28])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2)
    for i in range(6):
        cell = table[(0, i)]
        cell.set_facecolor('#3c8dbc')
        cell.set_text_props(weight='bold', color='white')
    ax.set_title('Table 2  Summary of thyroid cancer glycosylation biomarkers', fontsize=13, fontweight='bold', pad=20)
    plt.savefig(f'{FIG_DIR}/table2_biomarker_summary_en.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("[OK] Table 2 generated")

if __name__ == "__main__":
    print("=" * 60)
    print("Generating English-labeled figures for SCI manuscript")
    print("=" * 60)
    fig1(); fig2(); fig3(); fig4(); fig5(); table1(); table2()
    print("=" * 60)
    print("All figures generated!")
