#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
从清理后的数据库生成Cloudflare Worker的内联JSON数据
"""

import sqlite3
import json
import os

def generate_worker_data():
    conn = sqlite3.connect('output/acupoint_gene.db')
    cursor = conn.cursor()
    
    data = {}
    
    # 1. Stats
    cursor.execute('SELECT COUNT(*) FROM articles')
    articles = cursor.fetchone()[0]
    cursor.execute("SELECT COUNT(*) FROM genes WHERE quality_flag IN ('valid', 'mapped_to_official')")
    valid_genes = cursor.fetchone()[0]
    cursor.execute("SELECT COUNT(*) FROM genes WHERE quality_flag = 'gene_family'")
    gene_families = cursor.fetchone()[0]
    cursor.execute('SELECT COUNT(*) FROM acupoints')
    acupoints = cursor.fetchone()[0]
    cursor.execute("SELECT COUNT(*) FROM diseases WHERE entity_category != 'merged_synonym'")
    diseases = cursor.fetchone()[0]
    cursor.execute('SELECT COUNT(*) FROM articles WHERE has_fulltext = 1')
    fulltext = cursor.fetchone()[0]
    cursor.execute('SELECT COUNT(*) FROM article_gene ag JOIN genes g ON ag.gene_id = g.id WHERE g.quality_flag IN ("valid", "mapped_to_official")')
    article_gene = cursor.fetchone()[0]
    cursor.execute('SELECT COUNT(*) FROM article_acupoint')
    article_acupoint = cursor.fetchone()[0]
    cursor.execute('SELECT COUNT(*) FROM article_disease')
    article_disease = cursor.fetchone()[0]
    
    data['stats'] = {
        'articles': articles,
        'genes': valid_genes,
        'gene_families': gene_families,
        'acupoints': acupoints,
        'diseases': diseases,
        'fulltext': fulltext,
        'article_gene': article_gene,
        'article_acupoint': article_acupoint,
        'article_disease': article_disease,
    }
    
    # 2. Year data
    cursor.execute('SELECT year, COUNT(*) FROM articles WHERE year IS NOT NULL GROUP BY year ORDER BY year')
    data['year_data'] = [{'year': r[0], 'count': r[1]} for r in cursor.fetchall()]
    
    # 3. Top genes (valid only)
    cursor.execute('''
        SELECT g.id, g.gene_symbol, COUNT(DISTINCT ag.article_id) as cnt
        FROM genes g JOIN article_gene ag ON g.id = ag.gene_id
        WHERE g.quality_flag IN ('valid', 'mapped_to_official')
        GROUP BY g.id ORDER BY cnt DESC LIMIT 50
    ''')
    data['top_genes'] = [{'id': r[0], 'gene': r[1], 'count': r[2]} for r in cursor.fetchall()]
    
    # 4. Top gene families
    cursor.execute('''
        SELECT g.id, g.gene_symbol, COUNT(DISTINCT ag.article_id) as cnt
        FROM genes g JOIN article_gene ag ON g.id = ag.gene_id
        WHERE g.quality_flag = 'gene_family'
        GROUP BY g.id ORDER BY cnt DESC
    ''')
    data['top_gene_families'] = [{'id': r[0], 'family': r[1], 'count': r[2]} for r in cursor.fetchall()]
    
    # 5. Top acupoints
    cursor.execute('''
        SELECT a.id, a.name_cn, a.code, a.meridian, COUNT(DISTINCT aa.article_id) as cnt
        FROM acupoints a LEFT JOIN article_acupoint aa ON a.id = aa.acupoint_id
        GROUP BY a.id ORDER BY cnt DESC LIMIT 50
    ''')
    data['top_acupoints'] = [{'id': r[0], 'name_cn': r[1], 'code': r[2], 'meridian': r[3], 'count': r[4]} for r in cursor.fetchall()]
    
    # 6. Top diseases
    cursor.execute('''
        SELECT d.id, d.name, d.entity_category, COUNT(DISTINCT ad.article_id) as cnt
        FROM diseases d JOIN article_disease ad ON d.id = ad.disease_id
        WHERE d.entity_category != 'merged_synonym'
        GROUP BY d.id ORDER BY cnt DESC LIMIT 30
    ''')
    data['top_diseases'] = [{'id': r[0], 'name': r[1], 'category': r[2], 'count': r[3]} for r in cursor.fetchall()]
    
    # 7. Acupoint list
    cursor.execute('SELECT id, name_cn, code, meridian FROM acupoints ORDER BY name_cn')
    data['acupoint_list'] = [{'id': r[0], 'name': r[1], 'code': r[2], 'meridian': r[3]} for r in cursor.fetchall()]
    
    # 8. Gene list (valid + family)
    cursor.execute('''
        SELECT id, gene_symbol, quality_flag, entity_type 
        FROM genes 
        WHERE quality_flag IN ('valid', 'mapped_to_official', 'gene_family')
        ORDER BY gene_symbol
    ''')
    data['gene_list'] = [{'id': r[0], 'gene': r[1], 'quality_flag': r[2], 'entity_type': r[3]} for r in cursor.fetchall()]
    
    # 9. Disease list
    cursor.execute('''
        SELECT id, name, entity_category FROM diseases 
        WHERE entity_category != 'merged_synonym'
        ORDER BY name
    ''')
    data['disease_list'] = [{'id': r[0], 'name': r[1], 'category': r[2]} for r in cursor.fetchall()]
    
    # 10. Articles
    cursor.execute('''
        SELECT id, pmid, title, journal, year, doi, has_fulltext, article_type, language
        FROM articles ORDER BY year DESC, id DESC
    ''')
    data['articles'] = []
    pmid_to_article = {}
    for r in cursor.fetchall():
        art = {'id': r[0], 'pmid': r[1], 'title': r[2], 'journal': r[3], 'year': r[4], 
               'doi': r[5], 'has_fulltext': r[6], 'article_type': r[7], 'language': r[8]}
        data['articles'].append(art)
        pmid_to_article[str(r[1])] = art
    data['pmid_to_article'] = pmid_to_article
    
    # 11. PMID mappings
    cursor.execute('''
        SELECT a.pmid, g.gene_symbol 
        FROM articles a
        JOIN article_gene ag ON a.id = ag.article_id
        JOIN genes g ON ag.gene_id = g.id
        WHERE g.quality_flag IN ('valid', 'mapped_to_official', 'gene_family')
    ''')
    pmid_genes = {}
    for pmid, gene in cursor.fetchall():
        pmid_genes.setdefault(str(pmid), []).append(gene)
    data['pmid_genes'] = pmid_genes
    
    cursor.execute('''
        SELECT a.pmid, ac.name_cn 
        FROM articles a
        JOIN article_acupoint aa ON a.id = aa.article_id
        JOIN acupoints ac ON aa.acupoint_id = ac.id
    ''')
    pmid_acupoints = {}
    for pmid, apt in cursor.fetchall():
        pmid_acupoints.setdefault(str(pmid), []).append({'name': apt})
    data['pmid_acupoints'] = pmid_acupoints
    
    cursor.execute('''
        SELECT a.pmid, d.name 
        FROM articles a
        JOIN article_disease ad ON a.id = ad.article_id
        JOIN diseases d ON ad.disease_id = d.id
        WHERE d.entity_category != 'merged_synonym'
    ''')
    pmid_diseases = {}
    for pmid, disease in cursor.fetchall():
        pmid_diseases.setdefault(str(pmid), []).append(disease)
    data['pmid_diseases'] = pmid_diseases
    
    # 12. Reverse mappings
    acupoint_pmids = {}
    gene_pmids = {}
    disease_pmids = {}
    
    for pmid, apts in pmid_acupoints.items():
        for apt in apts:
            acupoint_pmids.setdefault(apt['name'], set()).add(int(pmid))
    for pmid, genes in pmid_genes.items():
        for gene in genes:
            gene_pmids.setdefault(gene, set()).add(int(pmid))
    for pmid, diseases in pmid_diseases.items():
        for disease in diseases:
            disease_pmids.setdefault(disease, set()).add(int(pmid))
    
    # Convert sets to lists for JSON serialization
    data['acupoint_pmids'] = {k: list(v) for k, v in acupoint_pmids.items()}
    data['gene_pmids'] = {k: list(v) for k, v in gene_pmids.items()}
    data['disease_pmids'] = {k: list(v) for k, v in disease_pmids.items()}
    
    # 13. Network data
    cursor.execute('''
        SELECT a.name_cn, g.gene_symbol, COUNT(*) as cnt
        FROM acupoints a
        JOIN article_acupoint aa ON a.id = aa.acupoint_id
        JOIN article_gene ag ON aa.article_id = ag.article_id
        JOIN genes g ON ag.gene_id = g.id
        WHERE g.quality_flag IN ('valid', 'mapped_to_official')
        GROUP BY a.id, g.id
        HAVING cnt >= 2
        ORDER BY cnt DESC
        LIMIT 300
    ''')
    data['network_data'] = [{'acupoint': r[0], 'gene': r[1], 'count': r[2]} for r in cursor.fetchall()]
    
    # 14. Journal stats
    cursor.execute('''
        SELECT journal, COUNT(*) as cnt FROM articles 
        WHERE journal IS NOT NULL AND journal != ''
        GROUP BY journal ORDER BY cnt DESC LIMIT 20
    ''')
    data['journal_stats'] = [{'journal': r[0], 'count': r[1]} for r in cursor.fetchall()]
    
    # 15. Gene functions (keep existing if available, otherwise empty)
    data['gene_functions'] = {}
    
    # 16. Acupoint info
    data['acupoint_info'] = {}
    for apt in data['acupoint_list']:
        data['acupoint_info'][apt['name']] = {
            'code': apt['code'],
            'meridian': apt['meridian']
        }
    
    conn.close()
    
    # Save as JSON file
    output_path = 'cloudflare-test-worker/db_data.json'
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(data, f, ensure_ascii=False)
    
    print(f"Worker data generated: {output_path}")
    print(f"  Articles: {data['stats']['articles']}")
    print(f"  Valid genes: {data['stats']['genes']}")
    print(f"  Gene families: {data['stats']['gene_families']}")
    print(f"  Acupoints: {data['stats']['acupoints']}")
    print(f"  Diseases: {data['stats']['diseases']}")
    print(f"  Network edges: {len(data['network_data'])}")
    
    return data


if __name__ == '__main__':
    generate_worker_data()
