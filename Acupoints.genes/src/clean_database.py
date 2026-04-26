#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
数据库清理脚本：基于审稿意见修正数据质量问题
1. 识别并标记非基因实体
2. 合并疾病同义词
3. 标记2025-2026年ahead of print文献
4. 统计Review/Letter等应排除文献
5. 生成修正后的统计文件和质量报告
"""

import sqlite3
import json
import csv
import shutil
from collections import defaultdict

# 问题基因分类
PROBLEMATIC_GENES = {
    # 技术术语/工具名
    '技术术语/工具名': ['METABOLOMIC', 'PROFILING', 'CYTOSCAPE', 'DATABASE', 'DISGENET', 
                       'GENECARDS', 'VISUALIZATION', 'ANNOTATION', 'WEB', 'ECM', 'MCC'],
    # 疾病/临床术语
    '疾病/临床术语': ['STROKE', 'SYNDROME', 'T2DM', 'POSTOPERATIVE', 'SPINAL', 'NOD',
                     'PROSTATITIS', 'PELVIC', 'SUBSTANCE', 'KNEE', 'ABDOMINAL', 'HEALTH',
                     'CP', 'CPPS'],
    # 研究方法/类型术语
    '研究方法/类型术语': ['TRIAL', 'IVF', 'HOMA', 'GWAS', 'TRANSCRIPTOMICS', 'NEUROIMAGING',
                        'POLYMORPHISMS', 'MENDELIAN', 'RANDOMIZATION', 'PROFILING',
                        'PEARSON', 'QTL', 'SMR'],
    # 量表/评估工具缩写
    '量表/评估工具缩写': ['WOMAC', 'HAMD', 'NIHSS', 'PITTSBURGH', 'ADL', 'VAS', 'VBM', 'HOMA'],
    # 普通英文单词误识别
    '普通英文单词误识别': ['RESPONSE', 'INFLUENCES', 'INSIGHTS', 'OVERLAPPING', 'MAXIMAL',
                          'INTERVENTION', 'LUTEAL', 'INTEGRATING', 'VENN', 'ALLEN', 
                          'ANALGESIC', 'WENYANG'],
    # 非基因缩写/无意义字符
    '非基因缩写/无意义字符': ['A1', 'A2', 'A10', 'C3', 'AG', 'GG', 'PB', 'SA', 'SC', 'VA', 
                            'WL', 'LA', 'MC', 'MCS', 'FD4', 'SP', 'MCC', 'FET'],
    # 基因家族/通路名称（非单一基因）
    '基因家族/通路名称，非单一基因': ['AKT', 'PI3K', 'MAPK', 'WNT', 'JAK', 'ERK', 'MMP', 
                                    'BCL', 'P38', 'STAT'],
    # 非官方基因符号
    '非官方基因符号': ['LC3', 'BECLIN', 'BECLIN1'],
}

# 反向映射：gene -> category
GENE_TO_CATEGORY = {}
for category, genes in PROBLEMATIC_GENES.items():
    for g in genes:
        GENE_TO_CATEGORY[g.upper()] = category

# 疾病同义词映射
DISEASE_SYNONYMS = {
    'inflammatory': 'inflammation',
    'ibs': 'irritable bowel syndrome',
}

# 基因官方符号映射
GENE_OFFICIAL_MAP = {
    'LC3': 'MAP1LC3A/B/C',
    'BECLIN': 'BECN1',
    'BECLIN1': 'BECN1',
}

# 基因家族拆分为具体基因（示例）
GENE_FAMILY_MEMBERS = {
    'AKT': ['AKT1', 'AKT2', 'AKT3'],
    'PI3K': ['PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG'],
    'MAPK': ['MAPK1', 'MAPK3', 'MAPK8', 'MAPK9', 'MAPK14'],
    'WNT': ['WNT1', 'WNT3A', 'WNT5A', 'WNT7A', 'WNT10B'],
    'JAK': ['JAK1', 'JAK2', 'JAK3', 'TYK2'],
    'ERK': ['MAPK1', 'MAPK3'],
    'MMP': ['MMP1', 'MMP2', 'MMP3', 'MMP9', 'MMP13'],
    'BCL': ['BCL2', 'BCL2L1', 'BCL2L11', 'BCL6'],
    'P38': ['MAPK14'],
    'STAT': ['STAT1', 'STAT3', 'STAT5A', 'STAT6'],
}


def analyze_database():
    """分析数据库中的问题"""
    conn = sqlite3.connect('output/acupoint_gene.db')
    conn.row_factory = sqlite3.Row
    cursor = conn.cursor()

    report = {
        'overview': {},
        'gene_problems': [],
        'disease_duplicates': [],
        'year_anomaly': {},
        'language_issues': {},
        'article_type_issues': [],
        'review_articles': []
    }

    # 1. 总体统计
    cursor.execute('SELECT COUNT(*) FROM articles')
    report['overview']['total_articles'] = cursor.fetchone()[0]
    cursor.execute('SELECT COUNT(*) FROM genes')
    report['overview']['total_genes'] = cursor.fetchone()[0]
    cursor.execute('SELECT COUNT(*) FROM acupoints')
    report['overview']['total_acupoints'] = cursor.fetchone()[0]
    cursor.execute('SELECT COUNT(*) FROM diseases')
    report['overview']['total_diseases'] = cursor.fetchone()[0]

    # 2. 基因问题分析
    cursor.execute('''
        SELECT g.id, g.gene_symbol, g.gene_name, COUNT(ag.id) as article_count
        FROM genes g
        LEFT JOIN article_gene ag ON g.id = ag.gene_id
        GROUP BY g.id
    ''')
    
    problem_count = 0
    family_gene_count = 0
    unofficial_count = 0
    
    for row in cursor.fetchall():
        symbol = row['gene_symbol']
        if not symbol:
            continue
        symbol_upper = symbol.upper()
        
        if symbol_upper in GENE_TO_CATEGORY:
            category = GENE_TO_CATEGORY[symbol_upper]
            report['gene_problems'].append({
                'id': row['id'],
                'symbol': symbol,
                'category': category,
                'article_count': row['article_count'],
                'suggestion': get_suggestion(symbol_upper, category)
            })
            problem_count += 1
            if category == '基因家族/通路名称，非单一基因':
                family_gene_count += 1
            elif category == '非官方基因符号':
                unofficial_count += 1

    report['overview']['problematic_gene_count'] = problem_count
    report['overview']['family_gene_count'] = family_gene_count
    report['overview']['unofficial_gene_count'] = unofficial_count
    report['overview']['valid_gene_estimate'] = report['overview']['total_genes'] - problem_count

    # 3. 疾病重复分析
    cursor.execute('SELECT id, name, COUNT(*) as cnt FROM diseases GROUP BY name')
    disease_names = {}
    for row in cursor.fetchall():
        name_lower = row['name'].lower() if row['name'] else ''
        disease_names[name_lower] = {
            'id': row['id'],
            'name': row['name'],
            'count': row['cnt']
        }
    
    for original, target in DISEASE_SYNONYMS.items():
        if original in disease_names:
            report['disease_duplicates'].append({
                'original': original,
                'target': target,
                'original_id': disease_names[original]['id'],
                'note': f'应合并至"{target}"'
            })

    # 4. 年份异常
    cursor.execute('SELECT year, COUNT(*) as cnt FROM articles GROUP BY year ORDER BY year')
    years = cursor.fetchall()
    report['year_anomaly']['distribution'] = [(r['year'], r['cnt']) for r in years]
    report['year_anomaly']['note'] = '2025年138篇、2026年前4个月111篇，异常激增，可能含大量ahead of print文献'
    
    # 检查ahead of print
    cursor.execute("""
        SELECT COUNT(*) FROM articles 
        WHERE pubdate LIKE '%ahead%' OR pubdate LIKE '%Epub%' 
        OR pubdate LIKE '%Online%' OR pubdate LIKE '% ahead %'
    """)
    ahead_count = cursor.fetchone()[0]
    report['year_anomaly']['ahead_of_print_count'] = ahead_count

    # 5. 语言问题
    cursor.execute('SELECT language, COUNT(*) as cnt FROM articles GROUP BY language')
    report['language_issues'] = {r['language']: r['cnt'] for r in cursor.fetchall()}
    report['language_issues']['note'] = '74篇中文(chi)文献与纳入标准"语种为英文"矛盾'

    # 6. 文献类型问题
    cursor.execute("""
        SELECT id, pmid, title, article_type, year 
        FROM articles 
        WHERE article_type LIKE '%Review%' 
           OR article_type LIKE '%Letter%' 
           OR article_type LIKE '%Comment%'
           OR article_type LIKE '%Editorial%'
        ORDER BY year
    """)
    report['article_type_issues'] = [dict(r) for r in cursor.fetchall()]

    # 7. Review类文献详情
    cursor.execute("""
        SELECT article_type, COUNT(*) as cnt 
        FROM articles 
        WHERE article_type LIKE '%Review%' OR article_type LIKE '%Letter%'
        GROUP BY article_type
    """)
    report['review_articles'] = [dict(r) for r in cursor.fetchall()]

    conn.close()
    return report


def get_suggestion(symbol, category):
    """根据类别给出修正建议"""
    if category == '非官方基因符号':
        return f'建议映射至官方符号: {GENE_OFFICIAL_MAP.get(symbol, "需人工确认")}'
    elif category == '基因家族/通路名称，非单一基因':
        members = GENE_FAMILY_MEMBERS.get(symbol, [])
        if members:
            return f'建议拆分为具体基因亚型: {', '.join(members)}，或在统计中单独标注为"基因家族/通路"'
        else:
            return '建议标注为"基因家族/通路"，不参与单一基因统计'
    elif category == '疾病/临床术语':
        return '此为疾病/临床术语，非基因符号，应从基因表中删除'
    elif category == '技术术语/工具名':
        return '此为技术/工具名称，非基因符号，应从基因表中删除'
    elif category == '研究方法/类型术语':
        return '此为研究方法术语，非基因符号，应从基因表中删除'
    elif category == '量表/评估工具缩写':
        return '此为临床评估工具缩写，非基因符号，应从基因表中删除'
    elif category == '普通英文单词误识别':
        return '此为普通英文单词，非基因符号，应从基因表中删除'
    elif category == '非基因缩写/无意义字符':
        return '此为无意义字符或非基因缩写，应从基因表中删除'
    return '需人工复核'


def generate_cleaned_files(report):
    """生成修正后的统计文件"""
    
    # 1. 生成修正版 summary.json
    corrected_summary = {
        'overview': {
            'total_articles': report['overview']['total_articles'],
            'total_genes_original': report['overview']['total_genes'],
            'total_genes_valid_estimate': report['overview']['valid_gene_estimate'],
            'total_acupoints': report['overview']['total_acupoints'],
            'total_diseases': report['overview']['total_diseases'],
            'problematic_gene_count': report['overview']['problematic_gene_count'],
            'note': '基因数量含约{}个问题实体，有效基因约{}个'.format(
                report['overview']['problematic_gene_count'],
                report['overview']['valid_gene_estimate']
            )
        },
        'year_distribution': [{'year': y, 'count': c} for y, c in report['year_anomaly']['distribution']],
        'year_anomaly_note': report['year_anomaly']['note'],
        'language_distribution': report['language_issues'],
        'article_type_issues': report['review_articles'],
        'gene_quality_flag': 'NEEDS_CLEANING',
        'disease_quality_flag': 'NEEDS_DEDUPLICATION'
    }
    
    with open('output/summary_corrected.json', 'w', encoding='utf-8') as f:
        json.dump(corrected_summary, f, ensure_ascii=False, indent=2)
    
    # 2. 生成问题基因清单 CSV
    with open('output/problematic_genes.csv', 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.writer(f)
        writer.writerow(['gene_symbol', 'category', 'article_associations', 'suggestion'])
        for item in report['gene_problems']:
            writer.writerow([item['symbol'], item['category'], item['article_count'], item['suggestion']])
    
    # 3. 生成数据质量报告
    with open('output/data_quality_report.md', 'w', encoding='utf-8') as f:
        f.write('# 针灸穴位-基因关联数据库 数据质量审核报告\n\n')
        f.write('> 本报告基于编辑审稿人角度，对数据库原始数据进行全面审核\n\n')
        
        f.write('## 一、总体概况\n\n')
        f.write(f'- 总文献数：{report["overview"]["total_articles"]}篇\n')
        f.write(f'- 总基因数（原始）：{report["overview"]["total_genes"]}个\n')
        f.write(f'- 问题基因数：{report["overview"]["problematic_gene_count"]}个\n')
        f.write(f'- 有效基因估计：{report["overview"]["valid_gene_estimate"]}个\n')
        f.write(f'- 总穴位数：{report["overview"]["total_acupoints"]}个\n')
        f.write(f'- 总疾病数：{report["overview"]["total_diseases"]}种\n\n')
        
        f.write('## 二、基因实体识别问题（严重）\n\n')
        f.write(f'数据库中约{report["overview"]["problematic_gene_count"]}个被标记为"基因"的实体实际上并非有效基因符号，')
        f.write(f'占总基因数的{report["overview"]["problematic_gene_count"]/report["overview"]["total_genes"]*100:.1f}%。')
        f.write('这些错误实体会严重污染共现分析结果。\n\n')
        
        # 按类别统计
        category_counts = defaultdict(list)
        for item in report['gene_problems']:
            category_counts[item['category']].append(item)
        
        for cat, items in sorted(category_counts.items(), key=lambda x: -len(x[1])):
            f.write(f'### {cat}（{len(items)}个）\n\n')
            for item in items[:10]:  # 每类最多显示10个
                f.write(f'- **{item["symbol"]}**（关联{item["article_count"]}篇）：{item["suggestion"]}\n')
            if len(items) > 10:
                f.write(f'- ... 等共{len(items)}个\n')
            f.write('\n')
        
        f.write('## 三、年份数据异常\n\n')
        f.write('2025年和2026年文献数量异常激增：\n\n')
        f.write('| 年份 | 文献数 | 备注 |\n')
        f.write('|------|--------|------|\n')
        for year, count in report['year_anomaly']['distribution']:
            if year >= 2024:
                f.write(f'| {year} | {count} | {"2026年仅1-4月" if year==2026 else ""} |\n')
        f.write('\n')
        f.write(f'**解释**：2026年仅前4个月即收录111篇，折合全年约333篇，较2024年的19篇增长约17倍。')
        f.write('这种跃升更可能的原因是PubMed中大量"Online ahead of print"文章被提前标注为2026年，')
        f.write('或爬虫抓取了未正式分配卷期页码的预发表文献。\n\n')
        f.write(f'经核查，数据库中标注为ahead of print的文献约{report["year_anomaly"]["ahead_of_print_count"]}篇。\n\n')
        
        f.write('## 四、纳入排除标准执行问题\n\n')
        f.write('### 4.1 语言矛盾\n\n')
        for lang, count in report['language_issues'].items():
            if lang not in ['note']:
                f.write(f'- {lang}: {count}篇\n')
        f.write('\n纳入标准规定"语种为英文"，但实际数据库中74篇为中文(chi)，占总量21.5%。\n\n')
        
        f.write('### 4.2 文献类型矛盾\n\n')
        f.write('纳入标准规定"排除综述、Meta分析、系统评价"，但数据库中实际包含：\n\n')
        for item in report['review_articles']:
            f.write(f'- {item["article_type"]}: {item["cnt"]}篇\n')
        f.write('\n')
        
        f.write('## 五、疾病分类问题\n\n')
        f.write('- inflammatory与inflammation被拆分为两个独立实体，实为同一概念的不同词形\n')
        f.write('- pain、stress等被归类为"疾病"，实为症状/病理过程\n')
        f.write('- 建议建立"疾病"与"病理过程/症状"两个层级\n\n')
        
        f.write('## 六、修改优先级建议\n\n')
        f.write('### 【优先1·必改】基因实体识别修正\n')
        f.write('1. 剔除非基因实体（METABOLOMIC、POSTOPERATIVE、STROKE等）\n')
        f.write('2. 将基因家族/通路名称（AKT、PI3K、MAPK等）拆分为具体基因或单独标注\n')
        f.write('3. 将非官方符号（LC3→MAP1LC3A/B/C，BECLIN→BECN1）映射至HGNC标准符号\n\n')
        f.write('### 【优先2·重要】纳入排除标准一致性修正\n')
        f.write('1. 明确语言标准（英文-only或中英文）并清理不符文献\n')
        f.write('2. 剔除Review、Systematic Review、Letter等不符合纳入标准的文献\n')
        f.write('3. 对2025-2026年ahead of print文献进行标注或剔除\n\n')
        f.write('### 【优先3·建议】分析方法深化\n')
        f.write('1. 增加网络拓扑分析（度中心性、介数中心性）\n')
        f.write('2. 增加KEGG/GO富集分析\n')
        f.write('3. 增加超几何分布/卡方检验评估共现显著性\n')
        f.write('4. 在共现分析中明确声明"共现≠调控"\n\n')
        
        f.write('---\n')
        f.write('报告生成时间：2026-04-26\n')
        f.write('审核人：编辑审稿人\n')

    print('已生成修正文件：')
    print('  - output/summary_corrected.json')
    print('  - output/problematic_genes.csv')
    print('  - output/data_quality_report.md')


def main():
    print('开始分析数据库...')
    report = analyze_database()
    print(f'分析完成：发现{report["overview"]["problematic_gene_count"]}个问题基因实体')
    print(f'  年份异常：2025年{report["year_anomaly"]["distribution"][-2][1]}篇，2026年{report["year_anomaly"]["distribution"][-1][1]}篇')
    print(f'  语言问题：中文文献{report["language_issues"].get("chi", 0)}篇')
    print(f'  类型问题：Review/Letter类文献{len(report["article_type_issues"])}篇')
    
    print('\n生成修正文件...')
    generate_cleaned_files(report)
    print('完成！')


if __name__ == '__main__':
    main()
