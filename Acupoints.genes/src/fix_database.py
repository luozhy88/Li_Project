#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
数据库全面清理脚本
基于审稿意见修正所有数据质量问题
"""

import sqlite3
import json
import csv
import shutil
from collections import defaultdict
from datetime import datetime

# ========== 配置：问题基因清单 ==========

# 需要直接删除的非基因实体
NON_GENE_ENTITIES = {
    # 技术术语/工具名
    'METABOLOMIC', 'PROFILING', 'CYTOSCAPE', 'DATABASE', 'DISGENET', 
    'GENECARDS', 'VISUALIZATION', 'ANNOTATION', 'WEB', 'ECM', 'MCC',
    # 疾病/临床术语
    'STROKE', 'SYNDROME', 'T2DM', 'POSTOPERATIVE', 'SPINAL', 'NOD',
    'PROSTATITIS', 'PELVIC', 'SUBSTANCE', 'KNEE', 'ABDOMINAL', 'HEALTH',
    'CP', 'CPPS',
    # 研究方法/类型术语
    'TRIAL', 'IVF', 'HOMA', 'GWAS', 'TRANSCRIPTOMICS', 'NEUROIMAGING',
    'POLYMORPHISMS', 'MENDELIAN', 'RANDOMIZATION', 'PEARSON', 'QTL', 'SMR',
    # 量表/评估工具缩写
    'WOMAC', 'HAMD', 'NIHSS', 'PITTSBURGH', 'ADL', 'VAS', 'VBM',
    # 普通英文单词误识别
    'RESPONSE', 'INFLUENCES', 'INSIGHTS', 'OVERLAPPING', 'MAXIMAL',
    'INTERVENTION', 'LUTEAL', 'INTEGRATING', 'VENN', 'ALLEN', 'ANALGESIC', 'WENYANG',
    # 非基因缩写/无意义字符
    'A1', 'A2', 'A10', 'C3', 'AG', 'GG', 'PB', 'SA', 'SC', 'VA', 
    'WL', 'LA', 'MC', 'MCS', 'FD4', 'SP', 'FET',
}

# 基因家族/通路 → 需要特殊处理（从基因统计中移除，单独建表或标记）
GENE_FAMILIES = {
    'AKT': {'type': 'kinase_family', 'members': ['AKT1', 'AKT2', 'AKT3']},
    'PI3K': {'type': 'kinase_family', 'members': ['PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG']},
    'MAPK': {'type': 'kinase_family', 'members': ['MAPK1', 'MAPK3', 'MAPK8', 'MAPK9', 'MAPK14']},
    'WNT': {'type': 'signaling_family', 'members': ['WNT1', 'WNT3A', 'WNT5A', 'WNT7A', 'WNT10B']},
    'JAK': {'type': 'kinase_family', 'members': ['JAK1', 'JAK2', 'JAK3', 'TYK2']},
    'ERK': {'type': 'kinase_family', 'members': ['MAPK1', 'MAPK3']},
    'MMP': {'type': 'proteinase_family', 'members': ['MMP1', 'MMP2', 'MMP3', 'MMP9', 'MMP13']},
    'BCL': {'type': 'protein_family', 'members': ['BCL2', 'BCL2L1', 'BCL2L11', 'BCL6']},
    'P38': {'type': 'kinase', 'members': ['MAPK14']},
    'STAT': {'type': 'transcription_family', 'members': ['STAT1', 'STAT3', 'STAT5A', 'STAT6']},
}

# 非官方基因符号 → 官方映射
GENE_SYMBOL_MAP = {
    'LC3': 'MAP1LC3A/B/C',
    'BECLIN': 'BECN1',
    'BECLIN1': 'BECN1',
}

# 疾病同义词合并
DISEASE_MERGE_MAP = {
    'inflammatory': 'inflammation',
    'ibs': 'irritable bowel syndrome',
}


def connect_db():
    return sqlite3.connect('output/acupoint_gene.db')


def add_quality_columns(conn):
    """为表添加数据质量标记字段"""
    cursor = conn.cursor()
    
    # 为genes表添加quality_flag字段
    try:
        cursor.execute("ALTER TABLE genes ADD COLUMN quality_flag TEXT DEFAULT 'valid'")
        cursor.execute("ALTER TABLE genes ADD COLUMN original_symbol TEXT")
        cursor.execute("ALTER TABLE genes ADD COLUMN entity_type TEXT DEFAULT 'gene'")
    except sqlite3.OperationalError:
        pass  # 列已存在
    
    # 为articles表添加pub_status字段
    try:
        cursor.execute("ALTER TABLE articles ADD COLUMN pub_status TEXT DEFAULT 'unknown'")
        cursor.execute("ALTER TABLE articles ADD COLUMN data_quality_note TEXT")
    except sqlite3.OperationalError:
        pass
    
    # 为diseases表添加canonical_name字段
    try:
        cursor.execute("ALTER TABLE diseases ADD COLUMN canonical_name TEXT")
        cursor.execute("ALTER TABLE diseases ADD COLUMN entity_category TEXT DEFAULT 'disease'")
    except sqlite3.OperationalError:
        pass
    
    conn.commit()
    print("[1/8] 质量标记字段添加完成")


def clean_non_gene_entities(conn):
    """清理非基因实体"""
    cursor = conn.cursor()
    
    # 获取所有基因
    cursor.execute("SELECT id, gene_symbol FROM genes")
    genes = cursor.fetchall()
    
    removed_count = 0
    removed_associations = 0
    
    for gene_id, symbol in genes:
        if symbol and symbol.upper() in NON_GENE_ENTITIES:
            # 先删除关联
            cursor.execute("DELETE FROM article_gene WHERE gene_id = ?", (gene_id,))
            removed_associations += cursor.rowcount
            
            # 标记基因（不删除，保留记录但标记为无效）
            cursor.execute("""
                UPDATE genes 
                SET quality_flag = 'removed_non_gene_entity', 
                    entity_type = 'non_gene',
                    original_symbol = ?
                WHERE id = ?
            """, (symbol, gene_id))
            removed_count += 1
    
    conn.commit()
    print(f"[2/8] 非基因实体清理完成：标记{removed_count}个问题实体，删除{removed_associations}条关联")
    return removed_count, removed_associations


def handle_gene_families(conn):
    """处理基因家族/通路名称"""
    cursor = conn.cursor()
    
    handled_count = 0
    
    for family_name, info in GENE_FAMILIES.items():
        cursor.execute("SELECT id FROM genes WHERE gene_symbol = ?", (family_name,))
        row = cursor.fetchone()
        if row:
            gene_id = row[0]
            members_str = ','.join(info['members'])
            
            # 更新基因记录，标记为家族/通路
            cursor.execute("""
                UPDATE genes 
                SET quality_flag = 'gene_family', 
                    entity_type = ?,
                    gene_name = ?,
                    original_symbol = ?
                WHERE id = ?
            """, (info['type'], f'{family_name}家族({members_str})', family_name, gene_id))
            
            handled_count += 1
    
    conn.commit()
    print(f"[3/8] 基因家族处理完成：标记{handled_count}个家族/通路名称")
    return handled_count


def map_unofficial_symbols(conn):
    """映射非官方基因符号到官方符号"""
    cursor = conn.cursor()
    
    mapped_count = 0
    
    for unofficial, official in GENE_SYMBOL_MAP.items():
        cursor.execute("SELECT id FROM genes WHERE gene_symbol = ?", (unofficial,))
        row = cursor.fetchone()
        if not row:
            continue
            
        source_id = row[0]
        
        # 检查目标官方符号是否已存在
        cursor.execute("SELECT id FROM genes WHERE gene_symbol = ?", (official,))
        target_row = cursor.fetchone()
        
        if target_row:
            target_id = target_row[0]
            # 将所有article_gene关联从source重定向到target
            cursor.execute("""
                SELECT article_id FROM article_gene WHERE gene_id = ?
            """, (source_id,))
            articles = [r[0] for r in cursor.fetchall()]
            
            for article_id in articles:
                cursor.execute("""
                    SELECT 1 FROM article_gene 
                    WHERE article_id = ? AND gene_id = ?
                """, (article_id, target_id))
                if not cursor.fetchone():
                    cursor.execute("""
                        UPDATE article_gene 
                        SET gene_id = ? 
                        WHERE article_id = ? AND gene_id = ?
                    """, (target_id, article_id, source_id))
                else:
                    cursor.execute("""
                        DELETE FROM article_gene 
                        WHERE article_id = ? AND gene_id = ?
                    """, (article_id, source_id))
            
            # 标记原记录为已合并
            cursor.execute("""
                UPDATE genes 
                SET quality_flag = 'merged_to_official', 
                    entity_type = 'deprecated_symbol',
                    gene_name = ?,
                    original_symbol = ?
                WHERE id = ?
            """, (f'已合并至官方符号:{official}', unofficial, source_id))
        else:
            # 目标不存在，直接更新原记录
            cursor.execute("""
                UPDATE genes 
                SET quality_flag = 'mapped_to_official', 
                    original_symbol = gene_symbol,
                    gene_symbol = ?,
                    gene_name = ?
                WHERE id = ?
            """, (official, f'原符号:{unofficial}', source_id))
        
        mapped_count += 1
    
    conn.commit()
    print(f"[4/8] 非官方符号映射完成：映射{mapped_count}个符号")
    return mapped_count


def merge_disease_synonyms(conn):
    """合并疾病同义词"""
    cursor = conn.cursor()
    
    merged_count = 0
    
    for original, target in DISEASE_MERGE_MAP.items():
        # 查找原始实体和目标实体
        cursor.execute("SELECT id FROM diseases WHERE name = ?", (original,))
        orig_row = cursor.fetchone()
        
        cursor.execute("SELECT id FROM diseases WHERE name = ?", (target,))
        target_row = cursor.fetchone()
        
        if orig_row and target_row:
            orig_id = orig_row[0]
            target_id = target_row[0]
            
            # 将所有指向original的关联重定向到target
            cursor.execute("""
                SELECT article_id FROM article_disease 
                WHERE disease_id = ?
            """, (orig_id,))
            articles_to_move = [r[0] for r in cursor.fetchall()]
            
            for article_id in articles_to_move:
                # 检查是否已存在target关联
                cursor.execute("""
                    SELECT 1 FROM article_disease 
                    WHERE article_id = ? AND disease_id = ?
                """, (article_id, target_id))
                if not cursor.fetchone():
                    cursor.execute("""
                        UPDATE article_disease 
                        SET disease_id = ? 
                        WHERE article_id = ? AND disease_id = ?
                    """, (target_id, article_id, orig_id))
                else:
                    # 已存在，直接删除重复
                    cursor.execute("""
                        DELETE FROM article_disease 
                        WHERE article_id = ? AND disease_id = ?
                    """, (article_id, orig_id))
            
            # 标记原始疾病为已合并
            cursor.execute("""
                UPDATE diseases 
                SET canonical_name = ?, 
                    entity_category = 'merged_synonym'
                WHERE id = ?
            """, (target, orig_id))
            
            merged_count += 1
    
    conn.commit()
    print(f"[5/8] 疾病同义词合并完成：合并{merged_count}组同义词")
    return merged_count


def classify_diseases(conn):
    """对疾病/症状进行分类"""
    cursor = conn.cursor()
    
    # 定义病理过程/症状（非真正疾病）
    pathological_processes = {
        'pain', 'stress', 'inflammation', 'inflammatory pain', 'neuropathic pain',
        'chronic pain', 'fatigue', 'cognitive impairment', 'withdrawal',
        'cardiovascular', 'urinary', 'headache', 'nausea', 'fever',
        'constipation', 'diarrhea', 'dyspepsia', 'insomnia', 'anxiety',
        'depression', 'inflammatory'
    }
    
    cursor.execute("SELECT id, name FROM diseases")
    classified = 0
    
    for row in cursor.fetchall():
        disease_id, name = row
        if name and name.lower() in pathological_processes:
            cursor.execute("""
                UPDATE diseases 
                SET entity_category = 'pathological_process_or_symptom'
                WHERE id = ?
            """, (disease_id,))
            classified += 1
    
    conn.commit()
    print(f"[6/8] 疾病分类完成：标记{classified}个病理过程/症状")
    return classified


def mark_ahead_of_print(conn):
    """标记 ahead of print 文献"""
    cursor = conn.cursor()
    
    # 基于pubdate字段判断
    cursor.execute("""
        UPDATE articles 
        SET pub_status = 'ahead_of_print',
            data_quality_note = 'PubMed标注为Online ahead of print或Epub ahead of print'
        WHERE pubdate LIKE '%ahead%' 
           OR pubdate LIKE '%Epub%' 
           OR pubdate LIKE '%Online%'
    """)
    ahead_marked = cursor.rowcount
    
    # 标记2025-2026年为需核查
    cursor.execute("""
        UPDATE articles 
        SET data_quality_note = COALESCE(data_quality_note || '; ', '') || '2025-2026年文献激增，建议核查是否为ahead of print'
        WHERE year IN (2025, 2026) AND (pub_status IS NULL OR pub_status = 'unknown')
    """)
    
    # 标记Review/Letter类文献
    cursor.execute("""
        UPDATE articles 
        SET data_quality_note = COALESCE(data_quality_note || '; ', '') || 'Article type为Review/Letter/Commentary，与纳入标准矛盾'
        WHERE article_type LIKE '%Review%' 
           OR article_type LIKE '%Letter%' 
           OR article_type LIKE '%Comment%'
           OR article_type LIKE '%Editorial%'
    """)
    review_marked = cursor.rowcount
    
    # 标记中文文献
    cursor.execute("""
        UPDATE articles 
        SET data_quality_note = COALESCE(data_quality_note || '; ', '') || '语言为中文，与纳入标准"英文-only"矛盾'
        WHERE language = 'chi'
    """)
    chi_marked = cursor.rowcount
    
    conn.commit()
    print(f"[7/8] 文献质量标记完成：ahead_of_print={ahead_marked}, Review/Letter={review_marked}, 中文={chi_marked}")
    return ahead_marked, review_marked, chi_marked


def recalculate_statistics(conn):
    """重新计算所有统计数据"""
    cursor = conn.cursor()
    
    stats = {}
    
    # 总体统计（仅统计有效基因）
    cursor.execute("SELECT COUNT(*) FROM articles")
    stats['total_articles'] = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(*) FROM genes WHERE quality_flag = 'valid' OR quality_flag = 'mapped_to_official'")
    stats['valid_genes'] = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(*) FROM genes WHERE quality_flag = 'gene_family'")
    stats['gene_families'] = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(*) FROM genes WHERE quality_flag = 'removed_non_gene_entity'")
    stats['removed_non_genes'] = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(*) FROM acupoints")
    stats['total_acupoints'] = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(*) FROM diseases WHERE entity_category != 'merged_synonym'")
    stats['total_diseases'] = cursor.fetchone()[0]
    
    # 年份分布
    cursor.execute('''
        SELECT year, COUNT(*) FROM articles 
        WHERE year IS NOT NULL AND year BETWEEN 2010 AND 2026
        GROUP BY year ORDER BY year
    ''')
    stats['year_distribution'] = [{'year': r[0], 'count': r[1]} for r in cursor.fetchall()]
    
    # Top 15 有效基因
    cursor.execute('''
        SELECT g.gene_symbol, COUNT(DISTINCT ag.article_id) as cnt
        FROM genes g 
        JOIN article_gene ag ON g.id = ag.gene_id
        WHERE g.quality_flag IN ('valid', 'mapped_to_official')
        GROUP BY g.id 
        ORDER BY cnt DESC 
        LIMIT 15
    ''')
    stats['top_valid_genes'] = [{'gene': r[0], 'count': r[1]} for r in cursor.fetchall()]
    
    # Top 10 基因家族/通路
    cursor.execute('''
        SELECT g.gene_symbol, COUNT(DISTINCT ag.article_id) as cnt
        FROM genes g 
        JOIN article_gene ag ON g.id = ag.gene_id
        WHERE g.quality_flag = 'gene_family'
        GROUP BY g.id 
        ORDER BY cnt DESC 
        LIMIT 10
    ''')
    stats['top_gene_families'] = [{'family': r[0], 'count': r[1]} for r in cursor.fetchall()]
    
    # Top 10 穴位
    cursor.execute('''
        SELECT a.name_cn, a.code, COUNT(DISTINCT aa.article_id) as cnt
        FROM acupoints a 
        LEFT JOIN article_acupoint aa ON a.id = aa.acupoint_id
        GROUP BY a.id 
        ORDER BY cnt DESC 
        LIMIT 10
    ''')
    stats['top_acupoints'] = [{'name': r[0], 'code': r[1], 'count': r[2]} for r in cursor.fetchall()]
    
    # 疾病统计（分层）
    cursor.execute('''
        SELECT d.name, d.entity_category, COUNT(DISTINCT ad.article_id) as cnt
        FROM diseases d 
        JOIN article_disease ad ON d.id = ad.disease_id
        WHERE d.entity_category != 'merged_synonym'
        GROUP BY d.id 
        ORDER BY cnt DESC 
        LIMIT 15
    ''')
    stats['top_diseases'] = [{'name': r[0], 'category': r[1], 'count': r[2]} for r in cursor.fetchall()]
    
    # 足三里相关有效基因
    cursor.execute('''
        SELECT g.gene_symbol, COUNT(*) as cnt
        FROM acupoints a
        JOIN article_acupoint aa ON a.id = aa.acupoint_id
        JOIN article_gene ag ON aa.article_id = ag.article_id
        JOIN genes g ON ag.gene_id = g.id
        WHERE a.name_cn = '足三里' AND g.quality_flag IN ('valid', 'mapped_to_official')
        GROUP BY g.id ORDER BY cnt DESC LIMIT 10
    ''')
    stats['st36_valid_genes'] = [{'gene': r[0], 'count': r[1]} for r in cursor.fetchall()]
    
    # 百会相关有效基因
    cursor.execute('''
        SELECT g.gene_symbol, COUNT(*) as cnt
        FROM acupoints a
        JOIN article_acupoint aa ON a.id = aa.acupoint_id
        JOIN article_gene ag ON aa.article_id = ag.article_id
        JOIN genes g ON ag.gene_id = g.id
        WHERE a.name_cn = '百会' AND g.quality_flag IN ('valid', 'mapped_to_official')
        GROUP BY g.id ORDER BY cnt DESC LIMIT 8
    ''')
    stats['gv20_valid_genes'] = [{'gene': r[0], 'count': r[1]} for r in cursor.fetchall()]
    
    # 共现对统计（基于有效基因）
    cursor.execute('''
        SELECT a.name_cn, g.gene_symbol, COUNT(*) as cnt
        FROM acupoints a
        JOIN article_acupoint aa ON a.id = aa.acupoint_id
        JOIN article_gene ag ON aa.article_id = ag.article_id
        JOIN genes g ON ag.gene_id = g.id
        WHERE g.quality_flag IN ('valid', 'mapped_to_official')
        GROUP BY a.id, g.id
        ORDER BY cnt DESC LIMIT 10
    ''')
    stats['top_acupoint_gene_pairs'] = [
        {'acupoint': r[0], 'gene': r[1], 'count': r[2]} for r in cursor.fetchall()
    ]
    
    # 质量标记统计
    cursor.execute('SELECT COUNT(*) FROM articles WHERE data_quality_note IS NOT NULL')
    stats['flagged_articles'] = cursor.fetchone()[0]
    
    cursor.execute("SELECT COUNT(*) FROM articles WHERE language = 'chi'")
    stats['chinese_articles'] = cursor.fetchone()[0]
    
    cursor.execute("""
        SELECT COUNT(*) FROM articles 
        WHERE article_type LIKE '%Review%' OR article_type LIKE '%Letter%'
    """)
    stats['review_letter_articles'] = cursor.fetchone()[0]
    
    return stats


def export_csv_files(conn, stats):
    """导出修正后的CSV文件"""
    cursor = conn.cursor()
    
    # 1. genes.csv - 有效基因
    cursor.execute('''
        SELECT gene_symbol, 
               CASE WHEN quality_flag = 'mapped_to_official' 
                    THEN original_symbol || '→' || gene_symbol 
                    ELSE gene_symbol END as display_name,
               COUNT(ag.id) as article_count
        FROM genes g
        LEFT JOIN article_gene ag ON g.id = ag.gene_id
        WHERE g.quality_flag IN ('valid', 'mapped_to_official')
        GROUP BY g.id
        ORDER BY article_count DESC
    ''')
    with open('output/genes_cleaned.csv', 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.writer(f)
        writer.writerow(['gene_symbol', 'official_note', 'article_count'])
        for row in cursor.fetchall():
            writer.writerow(row)
    
    # 2. gene_families.csv
    cursor.execute('''
        SELECT gene_symbol, gene_name, COUNT(ag.id) as article_count
        FROM genes g
        LEFT JOIN article_gene ag ON g.id = ag.gene_id
        WHERE g.quality_flag = 'gene_family'
        GROUP BY g.id
        ORDER BY article_count DESC
    ''')
    with open('output/gene_families.csv', 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.writer(f)
        writer.writerow(['family_name', 'description', 'article_count'])
        for row in cursor.fetchall():
            writer.writerow(row)
    
    # 3. acupoints.csv
    cursor.execute('''
        SELECT a.name_cn, a.code, a.meridian, COUNT(DISTINCT aa.article_id) as cnt
        FROM acupoints a
        LEFT JOIN article_acupoint aa ON a.id = aa.acupoint_id
        GROUP BY a.id
        ORDER BY cnt DESC
    ''')
    with open('output/acupoints_cleaned.csv', 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.writer(f)
        writer.writerow(['name_cn', 'code', 'meridian', 'article_count'])
        for row in cursor.fetchall():
            writer.writerow(row)
    
    # 4. diseases.csv（分层）
    cursor.execute('''
        SELECT d.name, d.entity_category, COUNT(DISTINCT ad.article_id) as cnt
        FROM diseases d
        JOIN article_disease ad ON d.id = ad.disease_id
        WHERE d.entity_category != 'merged_synonym'
        GROUP BY d.id
        ORDER BY cnt DESC
    ''')
    with open('output/diseases_cleaned.csv', 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.writer(f)
        writer.writerow(['name', 'category', 'article_count'])
        for row in cursor.fetchall():
            writer.writerow(row)
    
    # 5. 质量标记文献清单
    cursor.execute('''
        SELECT pmid, title, year, language, article_type, data_quality_note
        FROM articles
        WHERE data_quality_note IS NOT NULL
        ORDER BY year
    ''')
    with open('output/flagged_articles.csv', 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.writer(f)
        writer.writerow(['pmid', 'title', 'year', 'language', 'article_type', 'quality_note'])
        for row in cursor.fetchall():
            writer.writerow(row)
    
    print("[8/8] CSV文件导出完成")


def generate_clean_summary(stats):
    """生成修正版 summary.json"""
    with open('output/summary_cleaned.json', 'w', encoding='utf-8') as f:
        json.dump(stats, f, ensure_ascii=False, indent=2)
    print("修正版 summary_cleaned.json 已生成")


def main():
    print("=" * 60)
    print("针灸穴位-基因关联数据库 全面清理脚本")
    print("=" * 60)
    print()
    
    conn = connect_db()
    
    # 1. 添加质量标记字段
    add_quality_columns(conn)
    
    # 2. 清理非基因实体
    removed_genes, removed_assoc = clean_non_gene_entities(conn)
    
    # 3. 处理基因家族
    family_count = handle_gene_families(conn)
    
    # 4. 映射非官方符号
    mapped_count = map_unofficial_symbols(conn)
    
    # 5. 合并疾病同义词
    merged_diseases = merge_disease_synonyms(conn)
    
    # 6. 疾病分类
    classified_diseases = classify_diseases(conn)
    
    # 7. 标记问题文献
    ahead_marked, review_marked, chi_marked = mark_ahead_of_print(conn)
    
    # 8. 重新计算统计
    print("\n重新计算统计数据...")
    stats = recalculate_statistics(conn)
    
    # 导出文件
    export_csv_files(conn, stats)
    generate_clean_summary(stats)
    
    conn.close()
    
    # 打印清理报告
    print("\n" + "=" * 60)
    print("数据库清理完成报告")
    print("=" * 60)
    print(f"总文献数：{stats['total_articles']}篇")
    print(f"有效基因数：{stats['valid_genes']}个（原声称2083个）")
    print(f"基因家族/通路数：{stats['gene_families']}个")
    print(f"已移除非基因实体：{stats['removed_non_genes']}个")
    print(f"总穴位数：{stats['total_acupoints']}个")
    print(f"总疾病/症状数：{stats['total_diseases']}种")
    print(f"被标记问题文献：{stats['flagged_articles']}篇")
    print(f"  - 中文文献：{stats['chinese_articles']}篇")
    print(f"  - Review/Letter类：{stats['review_letter_articles']}篇")
    print()
    print("生成的文件：")
    print("  - output/summary_cleaned.json")
    print("  - output/genes_cleaned.csv")
    print("  - output/gene_families.csv")
    print("  - output/acupoints_cleaned.csv")
    print("  - output/diseases_cleaned.csv")
    print("  - output/flagged_articles.csv")
    print("  - output/acupoint_gene.db（已修改，原备份为 .backup）")
    print()
    print("=" * 60)


if __name__ == '__main__':
    main()
