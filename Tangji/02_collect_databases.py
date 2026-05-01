#!/usr/bin/env python3
"""
02_collect_databases.py
收集糖组学相关数据库和平台信息，输出到 input/databases/
"""

import json
import csv
import os
from datetime import datetime

OUTPUT_DIR = "input/databases"
os.makedirs(OUTPUT_DIR, exist_ok=True)

DATABASES = [
    {
        "id": "DB001",
        "name": "GlycoPOST",
        "url": "https://glycopost.glycosmos.org/",
        "type": "质谱原始数据存储库",
        "description": "糖组学质谱原始数据标准存储库，遵循FAIR原则和MIRAGE指南",
        "features": [
            "接受各种质谱仪的原始数据",
            "分配可追溯的登录号",
            "支持元数据输入（实验条件、仪器设置）",
            "与jPOST和GlyCosmos门户集成"
        ],
        "thyroid_cancer_data": "有甲状腺癌相关数据，如GPST000197（髓样甲状腺癌糖鞘脂）、GPST000495（PTC N-聚糖）",
        "access": "免费",
        "last_update": "持续更新",
        "organization": "GlyCosmos Project / 日本",
        "citation": "Watanabe Y et al. Nucleic Acids Res. 2021"
    },
    {
        "id": "DB002",
        "name": "UniCarbKB",
        "url": "https://www.unicarbkb.org/",
        "type": "糖蛋白知识库",
        "description": "国际合作的糖蛋白糖链结构数据库，整合了GlycoSuiteDB和EUROCarbDB",
        "features": [
            "糖链结构详细信息（组成、连接方式、序列）",
            "糖蛋白全局和位点特异性附着信息",
            "与UniProtKB和GlycoMod集成",
            "支持通过结构基序、单糖组成搜索"
        ],
        "thyroid_cancer_data": "包含甲状腺球蛋白(Tg)等甲状腺相关糖蛋白信息，但疾病特异性整合度低",
        "access": "免费",
        "last_update": "持续更新",
        "organization": "国际合作项目",
        "citation": "Campbell MP et al. Nucleic Acids Res. 2014"
    },
    {
        "id": "DB003",
        "name": "GlycoStore",
        "url": "https://www.glycostore.org/",
        "type": "色谱和质谱数据库",
        "description": " curated色谱、电泳和质谱组成数据库，从GlycoBase发展而来",
        "features": [
            "850种独特糖链结构（GlycoCT和SNFG格式）",
            "支持HILIC-U/HPLC、RP-U/HPLC、PGC色谱、CE-LIF数据",
            "可按实验类型、样本名称、相关糖蛋白、分类学搜索"
        ],
        "thyroid_cancer_data": "无特异性甲状腺癌数据模块",
        "access": "免费",
        "last_update": "2019年",
        "organization": "英国/爱尔兰",
        "citation": "Staples GO et al."
    },
    {
        "id": "DB004",
        "name": "GlyTouCan",
        "url": "https://glytoucan.org/",
        "type": "国际糖链结构库",
        "description": "国际糖链结构注册中心，为单个糖链分配登录号",
        "features": [
            "唯一标识符系统",
            "与GlycoPOST和UniCarb-DR关联",
            "支持糖链结构标准化表示"
        ],
        "thyroid_cancer_data": "通用数据库，无疾病特异性",
        "access": "免费",
        "last_update": "持续更新",
        "organization": "国际糖组学组织",
        "citation": "Aoki-Kinoshita KF et al."
    },
    {
        "id": "DB005",
        "name": "TCGA (The Cancer Genome Atlas)",
        "url": "https://portal.gdc.cancer.gov/",
        "type": "癌症多组学数据库",
        "description": "美国NCI的癌症基因组图谱项目，包含甲状腺癌(THCA)的基因组、转录组数据",
        "features": [
            "mRNA转录组数据（509例肿瘤，58例正常）",
            "临床特征信息",
            "体细胞突变数据",
            "生存和预后数据"
        ],
        "thyroid_cancer_data": "567例甲状腺癌样本，但糖基化数据整合度极低，主要为基因组数据",
        "access": "免费",
        "last_update": "持续更新",
        "organization": "NCI/NIH (美国)",
        "citation": "TCGA Research Network"
    },
    {
        "id": "DB006",
        "name": "THPA (The Human Protein Atlas)",
        "url": "https://www.proteinatlas.org/",
        "type": "蛋白质组学数据库",
        "description": "人类蛋白质图谱，包含正常和肿瘤组织的免疫组化图像",
        "features": [
            "正常甲状腺和甲状腺癌组织的IHC图像",
            "糖基转移酶表达数据",
            "细胞、组织和病理学信息"
        ],
        "thyroid_cancer_data": "包含HSPA5、KIF20A、SDC2等基因在甲状腺癌中的表达，但缺乏糖基化特异性数据",
        "access": "免费",
        "last_update": "持续更新",
        "organization": "瑞典",
        "citation": "Uhlen M et al. Science. 2015"
    },
    {
        "id": "DB007",
        "name": "GlyConnect",
        "url": "https://glyconnect.expasy.org/",
        "type": "糖蛋白组学平台",
        "description": "整合多种信息源帮助表征蛋白质糖基化分子组分的平台",
        "features": [
            "糖蛋白信息整合",
            "糖链结构注释",
            "与ExPASy集成"
        ],
        "thyroid_cancer_data": "无特异性甲状腺癌模块",
        "access": "免费",
        "last_update": "持续更新",
        "organization": "SIB (瑞士)",
        "citation": ""
    },
    {
        "id": "DB008",
        "name": "GlycoProtDB (GPDB)",
        "url": "https://acgg.asia/gpdb2",
        "type": "糖蛋白数据库",
        "description": "提供 curated 实验支持的糖蛋白信息，涵盖线虫、人类和小鼠",
        "features": [
            "糖基化位点信息",
            "糖肽获取过程",
            "LC/MS基于的糖蛋白组学技术"
        ],
        "thyroid_cancer_data": "人类数据可用，但无疾病特异性分类",
        "access": "免费",
        "last_update": "2019年5月",
        "organization": "ACGG Asia",
        "citation": ""
    },
    {
        "id": "DB009",
        "name": "UniPep / N-GlycositeAtlas",
        "url": "http://www.unipep.org/",
        "type": "N-糖基化位点数据库",
        "description": "人类N-连接糖基化位点数据库，聚焦于可作为生物标志物的糖肽",
        "features": [
            "1552个独特N-糖基化位点",
            "不同体液来源",
            "肼化学法捕获N-糖蛋白"
        ],
        "thyroid_cancer_data": "通用人类数据，无甲状腺癌特异性",
        "access": "免费",
        "last_update": "",
        "organization": "",
        "citation": ""
    },
    {
        "id": "DB010",
        "name": "CSDB (Carbohydrate Structure Database)",
        "url": "http://csdb.glycoscience.ru/database/",
        "type": "碳水化合物结构数据库",
        "description": "原核生物、植物和真菌的天然碳水化合物及其衍生物的结构和NMR信号数据库",
        "features": [
            "25397种化合物",
            "13151种生物体",
            "SFNG和BCSDB格式",
            "NMR光谱解释和预测工具"
        ],
        "thyroid_cancer_data": "无特异性数据",
        "access": "免费",
        "last_update": "2021年1月",
        "organization": "俄罗斯科学院",
        "citation": "Toukach PV et al."
    }
]

def save_json():
    path = os.path.join(OUTPUT_DIR, "glycosylation_databases.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump({
            "metadata": {
                "title": "糖组学相关数据库和平台收集",
                "created": datetime.now().isoformat(),
                "count": len(DATABASES),
                "note": "目前尚无专门针对甲状腺癌糖基化的交互分析平台"
            },
            "databases": DATABASES
        }, f, ensure_ascii=False, indent=2)
    print(f"[OK] JSON saved: {path}")

def save_csv():
    path = os.path.join(OUTPUT_DIR, "glycosylation_databases.csv")
    keys = ["id", "name", "url", "type", "description", "thyroid_cancer_data", "access", "last_update", "organization", "citation"]
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for item in DATABASES:
            row = {k: item.get(k, "") for k in keys}
            writer.writerow(row)
    print(f"[OK] CSV saved: {path}")

def save_markdown():
    path = os.path.join(OUTPUT_DIR, "glycosylation_databases.md")
    lines = [
        "# 糖组学相关数据库和平台汇总",
        "",
        f"> 收集时间: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"> 数据库总数: {len(DATABASES)} 个",
        f"> **重要发现**: 目前尚无专门针对甲状腺癌糖基化的交互分析平台！",
        "",
        "---",
        ""
    ]
    for item in DATABASES:
        lines.append(f"## {item['id']}. {item['name']}")
        lines.append("")
        lines.append(f"- **类型**: {item['type']}")
        lines.append(f"- **URL**: {item['url']}")
        lines.append(f"- **描述**: {item['description']}")
        lines.append(f"- **访问权限**: {item['access']}")
        lines.append(f"- **最近更新**: {item['last_update']}")
        lines.append(f"- **机构**: {item['organization']}")
        lines.append(f"- **引用**: {item['citation']}")
        lines.append(f"- **甲状腺癌数据**: {item['thyroid_cancer_data']}")
        lines.append("")
        lines.append("**主要功能**:")
        for feat in item['features']:
            lines.append(f"- {feat}")
        lines.append("")
        lines.append("---")
        lines.append("")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"[OK] Markdown saved: {path}")

if __name__ == "__main__":
    save_json()
    save_csv()
    save_markdown()
    print(f"\n数据库收集完成，共 {len(DATABASES)} 个。")
