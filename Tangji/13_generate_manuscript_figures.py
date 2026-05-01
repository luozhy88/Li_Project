#!/usr/bin/env python3
"""
13_generate_manuscript_figures.py
生成论文所需的图表，使用matplotlib和真实数据
"""

import sqlite3
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import os

plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

DB_PATH = "output/thyroid_glyco_db.sqlite"
FIG_DIR = "manuscript/figures"
TABLE_DIR = "manuscript/tables"
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TABLE_DIR, exist_ok=True)

con = sqlite3.connect(DB_PATH)

# ========================================================================
# 图1：甲状腺癌糖基化研究文献年份分布与趋势
# ========================================================================
def fig1_literature_trend():
    cursor = con.cursor()
    cursor.execute("SELECT year, COUNT(*) FROM studies WHERE year IS NOT NULL GROUP BY year ORDER BY year")
    rows = cursor.fetchall()
    years = [r[0] for r in rows]
    counts = [r[1] for r in rows]
    
    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(years, counts, color='#3c8dbc', edgecolor='white', width=0.6)
    
    # 添加数值标签
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 0.05,
                f'{int(count)}', ha='center', va='bottom', fontsize=10)
    
    ax.set_xlabel('发表年份', fontsize=12)
    ax.set_ylabel('文献数量（篇）', fontsize=12)
    ax.set_title('图1  甲状腺癌糖基化相关研究文献年份分布', fontsize=14, fontweight='bold')
    ax.set_xticks(years)
    ax.set_ylim(0, max(counts) + 1)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/fig1_literature_trend.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{FIG_DIR}/fig1_literature_trend.pdf', bbox_inches='tight')
    plt.close()
    print("[OK] 图1 生成完成")

# ========================================================================
# 图2：关键生物标志物诊断性能对比
# ========================================================================
def fig2_biomarker_auc():
    cursor = con.cursor()
    cursor.execute("SELECT name, performance_auc, performance_sensitivity, performance_specificity FROM biomarkers WHERE performance_auc IS NOT NULL")
    rows = cursor.fetchall()
    
    names = []
    aucs = []
    for r in rows:
        names.append(r[0][:30])
        aucs.append(r[1])
    
    fig, ax = plt.subplots(figsize=(8, 5))
    colors = ['#dd4b39' if a >= 0.9 else '#f39c12' if a >= 0.8 else '#999999' for a in aucs]
    bars = ax.barh(range(len(names)), aucs, color=colors, edgecolor='white', height=0.5)
    
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=9)
    ax.set_xlabel('AUC', fontsize=12)
    ax.set_title('图2  已报道甲状腺癌糖基化标志物诊断性能比较', fontsize=14, fontweight='bold')
    ax.axvline(x=0.8, color='orange', linestyle='--', linewidth=1, label='AUC=0.80')
    ax.axvline(x=0.9, color='green', linestyle='--', linewidth=1, label='AUC=0.90')
    ax.set_xlim(0.7, 1.0)
    
    for i, (bar, auc) in enumerate(zip(bars, aucs)):
        ax.text(auc + 0.01, i, f'{auc:.3f}', va='center', fontsize=9)
    
    ax.legend(loc='lower right', fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/fig2_biomarker_auc.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{FIG_DIR}/fig2_biomarker_auc.pdf', bbox_inches='tight')
    plt.close()
    print("[OK] 图2 生成完成")

# ========================================================================
# 图3：TCGA糖基转移酶差异表达（真实数据）
# ========================================================================
def fig3_tcga_glycogene_expression():
    cursor = con.cursor()
    cursor.execute("SELECT gene_symbol, normal_median_rpkm, ptc_median_rpkm, fold_change, significant FROM tcga_glycogene_expression ORDER BY fold_change DESC")
    rows = cursor.fetchall()
    
    genes = [r[0] for r in rows]
    normal = [r[1] for r in rows]
    ptc = [r[2] for r in rows]
    sig = [r[4] for r in rows]
    
    x = np.arange(len(genes))
    width = 0.35
    
    fig, ax = plt.subplots(figsize=(10, 6))
    bars1 = ax.bar(x - width/2, normal, width, label='Normal (n=20)', color='#00a65a', edgecolor='white')
    bars2 = ax.bar(x + width/2, ptc, width, label='PTC (n=20)', color='#dd4b39', edgecolor='white')
    
    # 显著性标记
    for i, s in enumerate(sig):
        if s == 1:
            max_h = max(normal[i], ptc[i])
            ax.plot([i-width/2, i+width/2], [max_h+0.3, max_h+0.3], 'k-', linewidth=1)
            ax.text(i, max_h+0.45, '*', ha='center', fontsize=14)
    
    ax.set_ylabel('Median RPKM', fontsize=12)
    ax.set_title('图3  TCGA数据库中甲状腺癌糖基转移酶表达差异（真实数据）\n数据来源：Bones et al. Cancers 2018', fontsize=13, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(genes, rotation=45, ha='right', fontsize=10)
    ax.legend(loc='upper right', fontsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/fig3_tcga_glycogene.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{FIG_DIR}/fig3_tcga_glycogene.pdf', bbox_inches='tight')
    plt.close()
    print("[OK] 图3 生成完成")

# ========================================================================
# 图4：平台架构设计图
# ========================================================================
def fig4_platform_architecture():
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 10)
    ax.axis('off')
    
    # 标题
    ax.text(6, 9.5, '图4  ThyGlycoPortal平台架构设计', fontsize=16, fontweight='bold', ha='center')
    
    # 颜色定义
    c_data = '#3498db'
    c_analysis = '#e74c3c'
    c_clinical = '#2ecc71'
    c_user = '#f39c12'
    
    def draw_box(x, y, w, h, text, color, fontsize=9):
        box = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.05", 
                              facecolor=color, edgecolor='black', alpha=0.8)
        ax.add_patch(box)
        ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=fontsize, 
                color='white', fontweight='bold', wrap=True)
    
    # 数据层
    draw_box(0.5, 7.5, 2.5, 1.2, '文献数据库\n(11项研究)', c_data)
    draw_box(3.5, 7.5, 2.5, 1.2, '糖链结构库\n(24种N-糖)', c_data)
    draw_box(6.5, 7.5, 2.5, 1.2, 'TCGA表达数据\n(20对配对)', c_data)
    draw_box(9.5, 7.5, 2.0, 1.2, '用户上传\n数据', c_data)
    
    ax.text(6, 7.0, '数据资源层', fontsize=11, ha='center', fontweight='bold')
    
    # 分析层
    draw_box(0.5, 4.8, 2.0, 1.2, '糖谱浏览\n与比较', c_analysis)
    draw_box(3.0, 4.8, 2.0, 1.2, '标志物筛选\n与验证', c_analysis)
    draw_box(5.5, 4.8, 2.0, 1.2, '糖基转移酶\n关联分析', c_analysis)
    draw_box(8.0, 4.8, 2.0, 1.2, '多组学\n整合', c_analysis)
    draw_box(10.2, 4.8, 1.5, 1.2, 'Meta\n分析', c_analysis)
    
    ax.text(6, 4.3, '核心分析层', fontsize=11, ha='center', fontweight='bold')
    
    # 临床决策层
    draw_box(2.0, 2.5, 3.0, 1.2, '列线图生成工具\n(LNM风险预测)', c_clinical)
    draw_box(6.0, 2.5, 3.0, 1.2, '临床决策支持\n(诊断/预后分层)', c_clinical)
    
    ax.text(6, 2.0, '临床决策层', fontsize=11, ha='center', fontweight='bold')
    
    # 用户层
    draw_box(4.0, 0.5, 4.0, 1.0, '终端用户界面 (R Shiny)', c_user, fontsize=11)
    
    # 箭头
    arrow_style = dict(arrowstyle='->', color='gray', lw=1.5)
    for x in [1.75, 4.75, 7.75, 10.5]:
        ax.annotate('', xy=(x, 6.0), xytext=(x, 7.5), arrowprops=arrow_style)
    
    for x in [1.5, 4.0, 6.5, 9.0, 10.95]:
        ax.annotate('', xy=(x, 3.7), xytext=(x, 4.8), arrowprops=arrow_style)
    
    ax.annotate('', xy=(6, 1.5), xytext=(6, 2.5), arrowprops=arrow_style)
    
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/fig4_platform_architecture.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{FIG_DIR}/fig4_platform_architecture.pdf', bbox_inches='tight')
    plt.close()
    print("[OK] 图4 生成完成")

# ========================================================================
# 图5：糖链变化模式示意图（基于文献）
# ========================================================================
def fig5_glycan_changes():
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    
    # Panel A: 血清糖组学复发预测
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
    ax.set_title('A. 复发预测\n(Kudelka 2023)', fontsize=11, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=9)
    ax.legend(fontsize=8)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Panel B: IgG BN诊断
    ax = axes[1]
    categories = ['HC\n(n=25)', 'TC\n(n=25)']
    bn = [18.5, 27.5]
    ax.bar(categories, bn, color=['#00a65a', '#dd4b39'], edgecolor='white', width=0.5)
    ax.set_ylabel('BN relative abundance (%)', fontsize=10)
    ax.set_title('B. IgG BN诊断\n(Zhang 2021, AUC=0.920)', fontsize=11, fontweight='bold')
    for i, v in enumerate(bn):
        ax.text(i, v + 0.5, f'{v}%', ha='center', fontsize=10, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Panel C: 糖基转移酶变化方向
    ax = axes[2]
    genes = ['ST6GAL1', 'MAN1A2', 'MAN2A1', 'B4GALT1', 'MGAT5']
    fc = [1.62, 0.69, 0.73, 0.74, 1.50]
    colors = ['#dd4b39' if f > 1 else '#00a65a' for f in fc]
    bars = ax.barh(range(len(genes)), fc, color=colors, edgecolor='white')
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=10)
    ax.set_xlabel('Fold change (PTC/Normal)', fontsize=10)
    ax.set_title('C. 糖基转移酶表达变化\n(Bones 2018, TCGA)', fontsize=11, fontweight='bold')
    ax.axvline(x=1, color='black', linestyle='--', linewidth=1)
    for i, (bar, f) in enumerate(zip(bars, fc)):
        ax.text(f + 0.05 if f > 1 else f - 0.05, i, f'{f:.2f}', 
                va='center', ha='left' if f > 1 else 'right', fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.suptitle('图5  甲状腺癌糖基化改变模式总结', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(f'{FIG_DIR}/fig5_glycan_changes.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{FIG_DIR}/fig5_glycan_changes.pdf', bbox_inches='tight')
    plt.close()
    print("[OK] 图5 生成完成")

# ========================================================================
# 表1：现有数据库对比表（输出为图片便于插入论文）
# ========================================================================
def table1_database_comparison():
    fig, ax = plt.subplots(figsize=(12, 5))
    ax.axis('tight')
    ax.axis('off')
    
    data = [
        ['数据库', '类型', '甲状腺癌特异性', '糖基化数据', '访问权限'],
        ['GlycoPOST', '质谱原始数据存储', '× 无专用模块', '✓ 原始质谱', '免费'],
        ['UniCarbKB', '糖蛋白知识库', '× 疾病整合度低', '✓ 糖链结构', '免费'],
        ['GlycoStore', '色谱/质谱数据库', '× 无', '✓ 实验数据', '免费'],
        ['TCGA', '癌症多组学', '△ 仅基因组', '× 极低', '免费'],
        ['THPA', '蛋白质图谱', '△ 仅IHC图像', '× 无', '免费'],
        ['GlyConnect', '糖蛋白组学平台', '× 无', '✓ 糖蛋白', '免费'],
        ['**ThyGlycoPortal(本研究)**', '**甲状腺癌糖基化专用**', '**✓ 专用平台**', '**✓ 整合分析**', '**免费**']
    ]
    
    table = ax.table(cellText=data, cellLoc='center', loc='center',
                     colWidths=[0.25, 0.25, 0.20, 0.15, 0.15])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2)
    
    # 表头样式
    for i in range(5):
        cell = table[(0, i)]
        cell.set_facecolor('#3c8dbc')
        cell.set_text_props(weight='bold', color='white')
    
    # 最后一行（本研究）高亮
    for i in range(5):
        cell = table[(7, i)]
        cell.set_facecolor('#d4edda')
        cell.set_text_props(weight='bold')
    
    ax.set_title('表1  现有糖组学/癌症数据库与ThyGlycoPortal比较', fontsize=13, fontweight='bold', pad=20)
    plt.savefig(f'{FIG_DIR}/table1_database_comparison.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{FIG_DIR}/table1_database_comparison.pdf', bbox_inches='tight')
    plt.close()
    print("[OK] 表1 生成完成")

# ========================================================================
# 表2：关键生物标志物汇总
# ========================================================================
def table2_biomarker_summary():
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.axis('tight')
    ax.axis('off')
    
    data = [
        ['标志物名称', '样本类型', '比较组', 'AUC', 'p值', '核心发现'],
        ['BN (bisecting neutral)', '血浆IgG', 'TC vs HC/BTN', '0.920', '<0.0001', '双分支型N-聚糖在TC中升高'],
        ['Glyco-panel', '血浆IgG', 'TC vs HC/BTN', '0.917', '<0.0001', '4种糖链组合分类'],
        ['G0F:G1F ratio', '血清', '复发 vs HC', '0.820', '0.004', '半乳糖基化降低预测复发'],
        ['CA4 + A2F0S0G', '血清', 'LNM vs NLNM', '0.702/0.658', '0.001/0.011', '多分支和半乳糖基化预测LNM'],
    ]
    
    table = ax.table(cellText=data, cellLoc='center', loc='center',
                     colWidths=[0.22, 0.12, 0.18, 0.10, 0.12, 0.26])
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 2)
    
    for i in range(6):
        cell = table[(0, i)]
        cell.set_facecolor('#3c8dbc')
        cell.set_text_props(weight='bold', color='white')
    
    ax.set_title('表2  甲状腺癌糖基化诊断/预后标志物性能汇总', fontsize=13, fontweight='bold', pad=20)
    plt.savefig(f'{FIG_DIR}/table2_biomarker_summary.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{FIG_DIR}/table2_biomarker_summary.pdf', bbox_inches='tight')
    plt.close()
    print("[OK] 表2 生成完成")

if __name__ == "__main__":
    print("=" * 60)
    print("生成论文图表")
    print("=" * 60)
    fig1_literature_trend()
    fig2_biomarker_auc()
    fig3_tcga_glycogene_expression()
    fig4_platform_architecture()
    fig5_glycan_changes()
    table1_database_comparison()
    table2_biomarker_summary()
    print("=" * 60)
    print("全部图表生成完成！")
    print(f"保存路径: {FIG_DIR}/")
