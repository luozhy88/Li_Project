#!/usr/bin/env python3
"""
03_collect_datasets.py
收集甲状腺癌糖基化相关可获取数据集，输出到 input/datasets/
"""

import json
import csv
import os
from datetime import datetime

OUTPUT_DIR = "input/datasets"
os.makedirs(OUTPUT_DIR, exist_ok=True)

DATASETS = [
    {
        "id": "DS001",
        "name": "GlycoPOST GPST000197",
        "repository": "GlycoPOST",
        "url": "https://glycopost.glycosmos.org/entry/GPST000197",
        "accession": "GPST000197",
        "description": "人甲状腺髓样癌糖鞘脂的LC-ESI/MS原始数据",
        "sample_info": {
            "type": "组织",
            "source": "Uppsala University Hospital, Sweden",
            "disease": "Medullary thyroid carcinoma (MTC)",
            "size": "有限数量"
        },
        "method": "LC-ESI/MS (LTQ linear quadrupole ion trap)",
        "data_format": "原始质谱文件 (.raw)",
        "availability": "公开",
        "mirage_compliant": "是",
        "citation": "Characterization of Human Medullary Thyroid Carcinoma Glycosphingolipids (Cancers 2021)"
    },
    {
        "id": "DS002",
        "name": "GlycoPOST GPST000495",
        "repository": "GlycoPOST",
        "url": "https://glycopost.glycosmos.org/entry/GPST000495",
        "accession": "GPST000495",
        "description": "甲状腺癌糖蛋白异常N-聚糖（含3-O-硫酸化半乳糖残基）的质谱数据",
        "sample_info": {
            "type": "组织",
            "source": "Beth Israel Deaconess Medical Center",
            "disease": "Papillary thyroid carcinoma (PTC)",
            "size": "配对PTC组织和正常组织"
        },
        "method": "HCD-MS/MS, MALDI-imaging",
        "data_format": "原始质谱文件",
        "availability": "公开",
        "mirage_compliant": "是",
        "citation": "Thyroid Carcinoma Glycoproteins Express Altered N-Glycans with 3-O-Sulfated Galactose Residues"
    },
    {
        "id": "DS003",
        "name": "TCGA-THCA RNAseq",
        "repository": "GDC Data Portal (TCGA)",
        "url": "https://portal.gdc.cancer.gov/",
        "accession": "TCGA-THCA",
        "description": "甲状腺癌的mRNA转录组数据，可用于糖基转移酶表达分析",
        "sample_info": {
            "type": "组织",
            "source": "多中心收集",
            "disease": "Thyroid carcinoma (THCA)",
            "size": "509例肿瘤组织，58例正常组织，共567例"
        },
        "method": "RNA sequencing",
        "data_format": "HTSeq-FPKM, counts",
        "availability": "公开",
        "mirage_compliant": "N/A",
        "citation": "TCGA Research Network",
        "note": "虽非直接糖组学数据，但可用于糖基转移酶(GALNT3, FUT8, ST6GAL1等)的表达关联分析"
    },
    {
        "id": "DS004",
        "name": "血浆IgG N-糖组诊断数据集 (PUMCH)",
        "repository": "文献原始数据",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8406750/",
        "accession": "N/A",
        "description": "北京协和医院收集的血浆IgG N-糖组数据，用于区分甲状腺癌与良性结节",
        "sample_info": {
            "type": "血浆",
            "source": "Peking Union Medical College Hospital",
            "disease": "Papillary thyroid carcinoma (PTC), Benign thyroid nodules (BTN)",
            "size": "发现队列：TC(n=25) vs HC(n=25)；验证队列：TC(n=47), BTN(n=18), HC(n=44)"
        },
        "method": "MALDI-TOF MS",
        "data_format": "文献中提供的相对定量数据",
        "availability": "需联系作者",
        "mirage_compliant": "未知",
        "citation": "Zhang ZJ et al. Front Oncol. 2021",
        "note": "IgG N-糖组诊断模型最成熟的标志物体系"
    },
    {
        "id": "DS005",
        "name": "血清N-糖组复发预测数据集",
        "repository": "文献原始数据",
        "url": "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10067050/",
        "accession": "N/A",
        "description": "复发性分化型甲状腺癌患者血清N-糖组MALDI-TOF MS数据",
        "sample_info": {
            "type": "血清",
            "source": "",
            "disease": "Recurrent differentiated thyroid cancer",
            "size": "13例复发患者 vs 15例健康对照"
        },
        "method": "MALDI-TOF MS (permethylation labeling)",
        "data_format": "文献中提供的35种N-聚糖组成和相对丰度",
        "availability": "需联系作者",
        "mirage_compliant": "未知",
        "citation": "Kudelka MR et al. Cancer Med. 2023",
        "note": "G0F:G1F比值作为复发预测标志物"
    },
    {
        "id": "DS006",
        "name": "PTMC血清N-糖组数据集",
        "repository": "文献原始数据",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9497917/",
        "accession": "N/A",
        "description": "乳头状甲状腺微小癌(PTMC)血清N-糖组数据，含LNM预测",
        "sample_info": {
            "type": "血清",
            "source": "",
            "disease": "Papillary thyroid microcarcinoma (PTMC) with/without LNM",
            "size": "多组比较：PTMC+LNM, PTMC-LNM, HC"
        },
        "method": "MALDI-TOF MS (Bruker rapifleXtreme)",
        "data_format": "96种直接检测N-糖特征和91种衍生特征",
        "availability": "需联系作者",
        "mirage_compliant": "未知",
        "citation": "Nomograms Based on Serum N-glycome for PTMC (Front Oncol. 2022)",
        "note": "包含诊断和LNM预测的列线图模型"
    },
    {
        "id": "DS007",
        "name": "血浆N-糖组BTN/TC/LNM数据集",
        "repository": "文献原始数据",
        "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8267918/",
        "accession": "N/A",
        "description": "血浆蛋白N-糖组数据，用于区分BTN和TC，并预测淋巴结转移",
        "sample_info": {
            "type": "血浆",
            "source": "",
            "disease": "Thyroid cancer (TC), Benign thyroid nodules (BTN)",
            "size": "TC, BTN, HC三组"
        },
        "method": "MALDI-TOF MS (含连接特异性唾液酸化信息)",
        "data_format": "文献中提供",
        "availability": "需联系作者",
        "mirage_compliant": "未知",
        "citation": "Front Endocrinol. 2021",
        "note": "包含连接特异性唾液酸化分析"
    }
]

def flatten_sample_info(ds):
    flat = dict(ds)
    si = flat.pop("sample_info", {})
    for k, v in si.items():
        flat[f"sample_{k}"] = v
    return flat

def save_json():
    path = os.path.join(OUTPUT_DIR, "thyroid_cancer_glycosylation_datasets.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump({
            "metadata": {
                "title": "甲状腺癌糖基化相关可获取数据集",
                "created": datetime.now().isoformat(),
                "count": len(DATASETS),
                "summary": {
                    "public_mass_spectrometry": 2,
                    "tcga_genomics": 1,
                    "literature_derived": 4
                }
            },
            "datasets": DATASETS
        }, f, ensure_ascii=False, indent=2)
    print(f"[OK] JSON saved: {path}")

def save_csv():
    path = os.path.join(OUTPUT_DIR, "thyroid_cancer_glycosylation_datasets.csv")
    flat_items = [flatten_sample_info(d) for d in DATASETS]
    keys = ["id", "name", "repository", "url", "accession", "description",
            "sample_type", "sample_source", "sample_disease", "sample_size",
            "method", "data_format", "availability", "mirage_compliant", "citation", "note"]
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for item in flat_items:
            row = {k: item.get(k, "") for k in keys}
            writer.writerow(row)
    print(f"[OK] CSV saved: {path}")

def save_markdown():
    path = os.path.join(OUTPUT_DIR, "thyroid_cancer_glycosylation_datasets.md")
    lines = [
        "# 甲状腺癌糖基化相关可获取数据集",
        "",
        f"> 收集时间: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"> 数据集总数: {len(DATASETS)} 个",
        f"> 公开质谱数据: 2 个 (GlycoPOST)",
        f"> 文献衍生数据: 4 个",
        f"> TCGA基因组数据: 1 个",
        "",
        "---",
        ""
    ]
    for item in DATASETS:
        lines.append(f"## {item['id']}. {item['name']}")
        lines.append("")
        lines.append(f"- **仓库**: {item['repository']}")
        lines.append(f"- **URL**: {item['url']}")
        lines.append(f"- **登录号**: {item.get('accession', 'N/A')}")
        lines.append(f"- **描述**: {item['description']}")
        lines.append(f"- **方法**: {item['method']}")
        lines.append(f"- **数据格式**: {item['data_format']}")
        lines.append(f"- **可用性**: {item['availability']}")
        lines.append(f"- **MIRAGE合规**: {item['mirage_compliant']}")
        lines.append(f"- **引用**: {item['citation']}")
        if item.get('note'):
            lines.append(f"- **备注**: {item['note']}")
        lines.append("")
        si = item.get('sample_info', {})
        lines.append("**样本信息**:")
        lines.append(f"- 类型: {si.get('type', 'N/A')}")
        lines.append(f"- 来源: {si.get('source', 'N/A')}")
        lines.append(f"- 疾病: {si.get('disease', 'N/A')}")
        lines.append(f"- 规模: {si.get('size', 'N/A')}")
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
    print(f"\n数据集收集完成，共 {len(DATASETS)} 个。")
