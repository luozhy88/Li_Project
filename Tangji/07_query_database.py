#!/usr/bin/env python3
"""
07_query_database.py
Phase 1: 甲状腺癌糖基化基础数据库查询与浏览工具
提供多种查询模式和数据导出功能
"""

import sqlite3
import json
import csv
import os
from tabulate import tabulate

DB_PATH = "output/thyroid_glyco_db.sqlite"
OUTPUT_DIR = "output/queries"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def get_conn():
    return sqlite3.connect(DB_PATH)

def query_all_studies():
    """查询所有研究概览"""
    conn = get_conn()
    cursor = conn.cursor()
    cursor.execute("""
        SELECT study_id, title, authors, year, method, sample_type, cancer_type, key_finding
        FROM studies ORDER BY year DESC
    """)
    rows = cursor.fetchall()
    conn.close()
    headers = ["ID", "Title", "Authors", "Year", "Method", "Sample", "Cancer Type", "Key Finding"]
    print("\n【所有研究概览】")
    print(tabulate(rows, headers=headers, tablefmt="grid", maxcolwidths=[6, 40, 15, 6, 15, 10, 20, 40]))
    return rows

def query_glycan_changes_by_study(study_id):
    """查询特定研究中的糖链变化趋势"""
    conn = get_conn()
    cursor = conn.cursor()
    cursor.execute("""
        SELECT 
            s.study_id,
            s.title,
            g.group_name,
            gs.snfg_name,
            gs.composition,
            gs.glycan_type,
            ga.abundance_value,
            ga.detection_method,
            s.sample_type
        FROM glycan_abundance ga
        JOIN samples s2 ON ga.sample_id = s2.sample_id
        JOIN clinical_groups g ON s2.group_id = g.group_id
        JOIN glycan_structures gs ON ga.glycan_id = gs.glycan_id
        JOIN studies s ON s2.study_id = s.study_id
        WHERE s.study_id = ?
        ORDER BY gs.composition
    """, (study_id,))
    rows = cursor.fetchall()
    conn.close()

    if not rows:
        print(f"\n[提示] 研究 {study_id} 暂无糖链丰度数据")
        return []

    headers = ["Study", "Group", "Glycan", "Composition", "Type", "Trend", "Method", "Sample"]
    print(f"\n【研究 {study_id} 的糖链变化趋势】")
    print(tabulate(rows, headers=headers, tablefmt="grid", maxcolwidths=[8, 15, 10, 12, 15, 8, 10, 10]))
    return rows

def query_biomarkers():
    """查询所有生物标志物"""
    conn = get_conn()
    cursor = conn.cursor()
    cursor.execute("""
        SELECT 
            b.biomarker_id,
            b.name,
            b.biomarker_type,
            b.sample_type,
            b.performance_auc,
            b.performance_sensitivity,
            b.performance_specificity,
            b.validation_status,
            s.title
        FROM biomarkers b
        JOIN studies s ON b.study_id = s.study_id
        ORDER BY b.performance_auc DESC
    """)
    rows = cursor.fetchall()
    conn.close()
    headers = ["ID", "Name", "Type", "Sample", "AUC", "Sens", "Spec", "Validation", "Study"]
    print("\n【生物标志物汇总】")
    print(tabulate(rows, headers=headers, tablefmt="grid", maxcolwidths=[6, 25, 12, 10, 6, 6, 6, 15, 30]))
    return rows

def query_glycan_structures(search_term=None):
    """查询糖链结构"""
    conn = get_conn()
    cursor = conn.cursor()
    if search_term:
        cursor.execute("""
            SELECT glycan_id, composition, snfg_name, mass, glycan_type, description
            FROM glycan_structures
            WHERE composition LIKE ? OR snfg_name LIKE ? OR description LIKE ?
            ORDER BY mass
        """, (f"%{search_term}%", f"%{search_term}%", f"%{search_term}%"))
    else:
        cursor.execute("""
            SELECT glycan_id, composition, snfg_name, mass, glycan_type, description
            FROM glycan_structures
            ORDER BY mass
        """)
    rows = cursor.fetchall()
    conn.close()
    headers = ["ID", "Composition", "SNFG Name", "Mass", "Type", "Description"]
    print(f"\n【糖链结构 {'(搜索: ' + search_term + ')' if search_term else ''}】")
    print(tabulate(rows, headers=headers, tablefmt="grid", maxcolwidths=[8, 12, 10, 8, 18, 35]))
    return rows

def query_enzyme_glycan_network():
    """查询糖基转移酶-糖链调控网络"""
    conn = get_conn()
    cursor = conn.cursor()
    cursor.execute("""
        SELECT 
            gt.gene_symbol,
            gt.enzyme_name,
            gt.enzyme_family,
            gt.pathway,
            gs.snfg_name,
            gs.composition,
            egl.linkage_type,
            egl.evidence_level
        FROM enzyme_glycan_links egl
        JOIN glycosyltransferases gt ON egl.enzyme_id = gt.enzyme_id
        JOIN glycan_structures gs ON egl.glycan_id = gs.glycan_id
        ORDER BY gt.gene_symbol
    """)
    rows = cursor.fetchall()
    conn.close()
    headers = ["Gene", "Enzyme", "Family", "Pathway", "Glycan", "Composition", "Linkage", "Evidence"]
    print("\n【糖基转移酶-糖链调控网络】")
    print(tabulate(rows, headers=headers, tablefmt="grid", maxcolwidths=[10, 25, 15, 15, 10, 12, 20, 10]))
    return rows

def query_cross_study_glycan_comparison(glycan_snfg):
    """跨研究比较特定糖链的变化趋势"""
    conn = get_conn()
    cursor = conn.cursor()
    cursor.execute("""
        SELECT 
            s.study_id,
            s.title,
            s.year,
            s.cancer_type,
            g.group_name,
            gs.snfg_name,
            gs.composition,
            ga.abundance_value,
            ga.detection_method
        FROM glycan_abundance ga
        JOIN samples s2 ON ga.sample_id = s2.sample_id
        JOIN clinical_groups g ON s2.group_id = g.group_id
        JOIN glycan_structures gs ON ga.glycan_id = gs.glycan_id
        JOIN studies s ON s2.study_id = s.study_id
        WHERE gs.snfg_name = ?
        ORDER BY s.year DESC, g.group_name
    """, (glycan_snfg,))
    rows = cursor.fetchall()
    conn.close()

    if not rows:
        print(f"\n[提示] 糖链 {glycan_snfg} 暂无跨研究数据")
        return []

    headers = ["Study", "Title", "Year", "Cancer", "Group", "Glycan", "Comp", "Trend", "Method"]
    print(f"\n【糖链 '{glycan_snfg}' 的跨研究比较】")
    print(tabulate(rows, headers=headers, tablefmt="grid", maxcolwidths=[8, 35, 6, 15, 15, 8, 10, 8, 10]))
    return rows

def query_database_stats():
    """数据库统计信息"""
    conn = get_conn()
    cursor = conn.cursor()
    stats = {}

    tables = ["studies", "glycan_structures", "clinical_groups", "samples",
              "glycan_abundance", "biomarkers", "glycosyltransferases", "enzyme_glycan_links"]
    for table in tables:
        cursor.execute(f"SELECT COUNT(*) FROM {table}")
        stats[table] = cursor.fetchone()[0]

    cursor.execute("SELECT DISTINCT cancer_type FROM studies")
    stats["cancer_types"] = [r[0] for r in cursor.fetchall()]

    cursor.execute("SELECT DISTINCT method FROM studies")
    stats["methods"] = [r[0] for r in cursor.fetchall()]

    cursor.execute("SELECT MIN(year), MAX(year) FROM studies")
    stats["year_range"] = cursor.fetchone()

    conn.close()
    return stats

def export_query_to_csv(query_name, rows, headers):
    """导出查询结果到CSV"""
    path = os.path.join(OUTPUT_DIR, f"{query_name}.csv")
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(rows)
    print(f"[OK] 导出 CSV: {path}")

def export_query_to_json(query_name, rows, headers):
    """导出查询结果到JSON"""
    path = os.path.join(OUTPUT_DIR, f"{query_name}.json")
    data = [dict(zip(headers, row)) for row in rows]
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
    print(f"[OK] 导出 JSON: {path}")

def interactive_menu():
    """交互式查询菜单"""
    while True:
        print("\n" + "=" * 50)
        print("  甲状腺癌糖基化数据库 - 查询工具")
        print("=" * 50)
        print("1. 查看所有研究概览")
        print("2. 查询特定研究的糖链变化")
        print("3. 查看所有生物标志物")
        print("4. 查询糖链结构")
        print("5. 查看酶-糖链调控网络")
        print("6. 跨研究比较特定糖链")
        print("7. 数据库统计信息")
        print("0. 退出")
        print("=" * 50)

        choice = input("\n请选择操作 [0-7]: ").strip()

        if choice == "1":
            query_all_studies()
        elif choice == "2":
            study_id = input("请输入研究ID (如 S001): ").strip()
            query_glycan_changes_by_study(study_id)
        elif choice == "3":
            query_biomarkers()
        elif choice == "4":
            term = input("请输入搜索关键词 (留空显示全部): ").strip()
            query_glycan_structures(term or None)
        elif choice == "5":
            query_enzyme_glycan_network()
        elif choice == "6":
            glycan = input("请输入糖链SNFG名称 (如 G0F, BN, CA4): ").strip()
            query_cross_study_glycan_comparison(glycan)
        elif choice == "7":
            stats = query_database_stats()
            print("\n【数据库统计信息】")
            for k, v in stats.items():
                print(f"  {k}: {v}")
        elif choice == "0":
            print("再见！")
            break
        else:
            print("无效选择，请重试。")

def auto_export_all():
    """自动导出所有核心查询结果"""
    print("\n【自动导出所有核心查询结果】")

    # 1. 所有研究
    conn = get_conn()
    cursor = conn.cursor()
    cursor.execute("SELECT * FROM studies ORDER BY year DESC")
    rows = cursor.fetchall()
    headers = [d[0] for d in cursor.description]
    export_query_to_csv("all_studies", rows, headers)
    export_query_to_json("all_studies", rows, headers)

    # 2. 所有糖链结构
    cursor.execute("SELECT * FROM glycan_structures ORDER BY mass")
    rows = cursor.fetchall()
    headers = [d[0] for d in cursor.description]
    export_query_to_csv("all_glycans", rows, headers)
    export_query_to_json("all_glycans", rows, headers)

    # 3. 所有生物标志物
    cursor.execute("""
        SELECT b.*, s.title as study_title FROM biomarkers b
        JOIN studies s ON b.study_id = s.study_id ORDER BY b.performance_auc DESC
    """)
    rows = cursor.fetchall()
    headers = [d[0] for d in cursor.description]
    export_query_to_csv("all_biomarkers", rows, headers)
    export_query_to_json("all_biomarkers", rows, headers)

    # 4. 糖链丰度趋势
    cursor.execute("""
        SELECT 
            s.study_id, s.title, g.group_name, gs.snfg_name, gs.composition,
            gs.glycan_type, ga.abundance_value, ga.detection_method
        FROM glycan_abundance ga
        JOIN samples s2 ON ga.sample_id = s2.sample_id
        JOIN clinical_groups g ON s2.group_id = g.group_id
        JOIN glycan_structures gs ON ga.glycan_id = gs.glycan_id
        JOIN studies s ON s2.study_id = s.study_id
        ORDER BY s.study_id, gs.composition
    """)
    rows = cursor.fetchall()
    headers = [d[0] for d in cursor.description]
    export_query_to_csv("glycan_abundance_trends", rows, headers)
    export_query_to_json("glycan_abundance_trends", rows, headers)

    # 5. 酶-糖链网络
    cursor.execute("""
        SELECT gt.gene_symbol, gt.enzyme_name, gt.enzyme_family, gt.pathway,
               gs.snfg_name, gs.composition, egl.linkage_type, egl.evidence_level
        FROM enzyme_glycan_links egl
        JOIN glycosyltransferases gt ON egl.enzyme_id = gt.enzyme_id
        JOIN glycan_structures gs ON egl.glycan_id = gs.glycan_id
        ORDER BY gt.gene_symbol
    """)
    rows = cursor.fetchall()
    headers = [d[0] for d in cursor.description]
    export_query_to_csv("enzyme_glycan_network", rows, headers)
    export_query_to_json("enzyme_glycan_network", rows, headers)

    conn.close()
    print("\n[OK] 所有核心查询结果已导出到 output/queries/")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "--export":
        auto_export_all()
    elif len(sys.argv) > 1 and sys.argv[1] == "--stats":
        stats = query_database_stats()
        print(json.dumps(stats, ensure_ascii=False, indent=2))
    else:
        interactive_menu()
