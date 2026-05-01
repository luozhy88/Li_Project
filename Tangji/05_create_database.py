#!/usr/bin/env python3
"""
05_create_database.py
Phase 1: 创建甲状腺癌糖基化基础数据库 (SQLite)
"""

import sqlite3
import os

DB_PATH = "output/thyroid_glyco_db.sqlite"
os.makedirs("output", exist_ok=True)

def create_database():
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    # 启用外键约束
    cursor.execute("PRAGMA foreign_keys = ON")

    # 1. 研究/文献表
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS studies (
        study_id TEXT PRIMARY KEY,
        title TEXT NOT NULL,
        authors TEXT,
        journal TEXT,
        year INTEGER,
        doi TEXT,
        pmid TEXT,
        pmcid TEXT,
        url TEXT,
        method TEXT,
        sample_type TEXT,
        cancer_type TEXT,
        key_finding TEXT,
        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
    )
    """)

    # 2. 糖链结构表
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS glycan_structures (
        glycan_id TEXT PRIMARY KEY,
        composition TEXT NOT NULL,
        snfg_name TEXT,
        mass REAL,
        glycan_type TEXT,
        description TEXT,
        glytoucan_id TEXT,
        structure_image_url TEXT
    )
    """)

    # 3. 临床分组表
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS clinical_groups (
        group_id TEXT PRIMARY KEY,
        study_id TEXT NOT NULL,
        group_name TEXT NOT NULL,
        group_description TEXT,
        sample_size INTEGER,
        age_mean REAL,
        age_range TEXT,
        gender_male_pct REAL,
        foreign key (study_id) references studies(study_id)
    )
    """)

    # 4. 样本表
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS samples (
        sample_id TEXT PRIMARY KEY,
        study_id TEXT NOT NULL,
        group_id TEXT NOT NULL,
        sample_type TEXT,
        source TEXT,
        notes TEXT,
        foreign key (study_id) references studies(study_id),
        foreign key (group_id) references clinical_groups(group_id)
    )
    """)

    # 5. 糖链丰度表 (核心数据表)
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS glycan_abundance (
        abundance_id INTEGER PRIMARY KEY AUTOINCREMENT,
        sample_id TEXT NOT NULL,
        glycan_id TEXT NOT NULL,
        abundance_value REAL,
        abundance_unit TEXT DEFAULT 'relative intensity',
        detection_method TEXT,
        foreign key (sample_id) references samples(sample_id),
        foreign key (glycan_id) references glycan_structures(glycan_id)
    )
    """)

    # 6. 生物标志物表
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS biomarkers (
        biomarker_id TEXT PRIMARY KEY,
        study_id TEXT NOT NULL,
        name TEXT NOT NULL,
        biomarker_type TEXT,
        description TEXT,
        sample_type TEXT,
        performance_auc REAL,
        performance_sensitivity REAL,
        performance_specificity REAL,
        cutoff_value TEXT,
        validation_status TEXT,
        foreign key (study_id) references studies(study_id)
    )
    """)

    # 7. 糖基转移酶表 (从TCGA等整合)
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS glycosyltransferases (
        enzyme_id TEXT PRIMARY KEY,
        gene_symbol TEXT NOT NULL,
        enzyme_name TEXT,
        enzyme_family TEXT,
        substrate_type TEXT,
        pathway TEXT,
        description TEXT
    )
    """)

    # 8. 酶-糖链关联表
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS enzyme_glycan_links (
        link_id INTEGER PRIMARY KEY AUTOINCREMENT,
        enzyme_id TEXT NOT NULL,
        glycan_id TEXT NOT NULL,
        linkage_type TEXT,
        evidence_level TEXT,
        foreign key (enzyme_id) references glycosyltransferases(enzyme_id),
        foreign key (glycan_id) references glycan_structures(glycan_id)
    )
    """)

    # 创建索引以提高查询效率
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_abundance_sample ON glycan_abundance(sample_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_abundance_glycan ON glycan_abundance(glycan_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_samples_group ON samples(group_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_groups_study ON clinical_groups(study_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_biomarkers_study ON biomarkers(study_id)")

    conn.commit()
    conn.close()
    print(f"[OK] 数据库创建成功: {DB_PATH}")
    print("[INFO] 包含8个表: studies, glycan_structures, clinical_groups, samples,")
    print("       glycan_abundance, biomarkers, glycosyltransferases, enzyme_glycan_links")

if __name__ == "__main__":
    create_database()
