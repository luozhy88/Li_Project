#!/usr/bin/env python3
"""
11_import_real_data_to_db.py
将文献提取的真实统计数据导入SQLite数据库
"""

import sqlite3
import json
import os

DB_PATH = "output/thyroid_glyco_db.sqlite"
JSON_PATH = "input/real_data/literature_extracted_real_data.json"

def get_conn():
    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA foreign_keys = ON")
    return conn

def create_real_data_tables():
    """创建真实数据存储表"""
    conn = get_conn()
    cursor = conn.cursor()
    
    # 文献报告的真实汇总统计表
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS literature_stats (
        stat_id INTEGER PRIMARY KEY AUTOINCREMENT,
        study TEXT NOT NULL,
        variable TEXT NOT NULL,
        group1 TEXT,
        group1_n INTEGER,
        group1_mean REAL,
        group1_sd REAL,
        group1_median REAL,
        group1_iqr_low REAL,
        group1_iqr_high REAL,
        group2 TEXT,
        group2_n INTEGER,
        group2_mean REAL,
        group2_sd REAL,
        group2_median REAL,
        group2_iqr_low REAL,
        group2_iqr_high REAL,
        p_value TEXT,
        direction TEXT,
        auc REAL,
        fold_change REAL,
        significant BOOLEAN,
        data_type TEXT DEFAULT 'mean_sd'
    )
    """)
    
    # TCGA糖基因真实表达数据 (Bones 2018)
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS tcga_glycogene_expression (
        gene_id INTEGER PRIMARY KEY AUTOINCREMENT,
        gene_symbol TEXT NOT NULL,
        enzyme_family TEXT,
        pathway TEXT,
        normal_median_rpkm REAL,
        ptc_median_rpkm REAL,
        fold_change REAL,
        p_value TEXT,
        direction TEXT,
        significant BOOLEAN,
        sample_size_normal INTEGER,
        sample_size_ptc INTEGER,
        citation TEXT
    )
    """)
    
    conn.commit()
    conn.close()
    print("[OK] 真实数据表创建完成")

def import_literature_stats():
    """导入文献汇总统计"""
    with open(JSON_PATH, "r", encoding="utf-8") as f:
        data = json.load(f)
    
    conn = get_conn()
    cursor = conn.cursor()
    
    records = []
    
    # Zhang 2021 - IgG glycan means
    z = data["datasets"]["zhang_2021_igg_glycan"]["discovery_cohort"]["glycans"]
    for k, v in z.items():
        records.append((
            "Zhang_2021", v["name"], "HC", 25, v.get("HC_mean"), v.get("HC_sd"),
            None, None, None,
            "TC", 25, v.get("TC_mean"), v.get("TC_sd"),
            None, None, None,
            v.get("p_value"), v.get("direction"), None, None, None, "mean_sd"
        ))
    
    # Zhang 2021 - BN biomarker
    bn = data["datasets"]["zhang_2021_igg_glycan"]["discovery_cohort"]["biomarkers"]["BN"]
    records.append((
        "Zhang_2021", bn["name"], "HC", 25, bn.get("HC_mean"), bn.get("HC_sd"),
        None, None, None,
        "TC", 25, bn.get("TC_mean"), bn.get("TC_sd"),
        None, None, None,
        bn.get("p_value"), "increased in TC", bn.get("AUC"), None, None, "mean_sd"
    ))
    
    # Kudelka 2023 - serum glycan means
    k = data["datasets"]["kudelka_2023_serum_recurrence"]["cohort"]["glycans_relative_abundance_percent"]
    for glycan_key, v in k.items():
        records.append((
            "Kudelka_2023", v["name"], "HC", 15, v.get("HC_mean"), v.get("HC_sd"),
            None, None, None,
            "Recurrent_DTC", 13, v.get("REC_mean"), v.get("REC_sd"),
            None, None, None,
            v.get("p_value"), v.get("direction"), None, None, None, "mean_sd"
        ))
    
    # Kudelka 2023 - G0F:G1F ratio (median/IQR)
    ratio = data["datasets"]["kudelka_2023_serum_recurrence"]["cohort"]["biomarker"]["G0F_G1F_ratio"]
    records.append((
        "Kudelka_2023", ratio["name"], "HC", 15, None, None,
        ratio.get("HC_median"), ratio.get("HC_IQR_low"), ratio.get("HC_IQR_high"),
        "Recurrent_DTC", 13, None, None,
        ratio.get("REC_median"), ratio.get("REC_IQR_low"), ratio.get("REC_IQR_high"),
        ratio.get("p_value"), "increased in recurrent", ratio.get("AUC"), None, None, "median_iqr"
    ))
    
    # PTMC 2022 - diagnostic markers (medians)
    p = data["datasets"]["ptmc_2022_nomogram"]["diagnostic_markers"]
    for k, v in p.items():
        records.append((
            "PTMC_2022", v["name"], "HC", None, None, None,
            v.get("HC_median"), None, None,
            "PTMC", None, None, None,
            v.get("PTMC_median"), None, None,
            v.get("p_value"), "altered in PTMC", None, None, None, "median"
        ))
    
    # PTMC 2022 - LNM prediction
    lnm = data["datasets"]["ptmc_2022_nomogram"]["LNM_prediction"]
    for k, v in lnm.items():
        records.append((
            "PTMC_2022_LNM", v.get("name", k), "NLNM", None, None, None,
            v.get("NLNM_median"), None, None,
            "LNM", None, None, None,
            v.get("LNM_median"), None, None,
            v.get("p_multivariate"), "associated with LNM", v.get("AUC"), None, None, "median"
        ))
    
    cursor.executemany("""
        INSERT INTO literature_stats
        (study, variable, group1, group1_n, group1_mean, group1_sd, group1_median, group1_iqr_low, group1_iqr_high,
         group2, group2_n, group2_mean, group2_sd, group2_median, group2_iqr_low, group2_iqr_high,
         p_value, direction, auc, fold_change, significant, data_type)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, records)
    
    conn.commit()
    conn.close()
    print(f"[OK] 导入 {len(records)} 条文献真实统计数据")

def import_tcga_glycogene_data():
    """导入TCGA糖基因真实表达数据"""
    with open(JSON_PATH, "r", encoding="utf-8") as f:
        data = json.load(f)
    
    conn = get_conn()
    cursor = conn.cursor()
    
    genes = data["datasets"]["bones_2018_tcga_glycogenes"]["glycogenes"]
    records = []
    
    for gene, info in genes.items():
        records.append((
            gene,
            info["family"],
            info["pathway"],
            info.get("normal_RPKM_median"),
            info.get("PTC_RPKM_median"),
            info.get("fold_change"),
            info.get("p_value"),
            info.get("direction"),
            1 if info.get("significant") else 0,
            20, 20,
            "Bones J et al. Cancers 2018 (TCGA RNAseq)"
        ))
    
    cursor.executemany("""
        INSERT INTO tcga_glycogene_expression
        (gene_symbol, enzyme_family, pathway, normal_median_rpkm, ptc_median_rpkm,
         fold_change, p_value, direction, significant, sample_size_normal, sample_size_ptc, citation)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, records)
    
    conn.commit()
    conn.close()
    print(f"[OK] 导入 {len(records)} 条TCGA糖基因真实表达数据")

def main():
    print("=" * 60)
    print("导入文献真实统计数据到数据库")
    print("=" * 60)
    create_real_data_tables()
    import_literature_stats()
    import_tcga_glycogene_data()
    print("=" * 60)
    print("完成！所有数据均为已发表文献的真实汇总统计量")
    print("=" * 60)

if __name__ == "__main__":
    main()
