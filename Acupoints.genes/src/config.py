"""
项目配置：针灸穴位-基因关联数据库
"""
import os

# 项目根目录
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# 输入/输出路径
INPUT_DIR = os.path.join(BASE_DIR, "input")
OUTPUT_DIR = os.path.join(BASE_DIR, "output")
LOGS_DIR = os.path.join(BASE_DIR, "logs")

# 子目录
PAPERS_DIR = os.path.join(INPUT_DIR, "papers")
XML_DIR = os.path.join(INPUT_DIR, "xml")
DB_PATH = os.path.join(OUTPUT_DIR, "acupoint_gene.db")

# NCBI API
NCBI_EMAIL = "research@example.com"  # 请替换为你的邮箱
NCBI_API_KEY = None  # 可选，申请API key可提高限速
NCBI_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
PMC_BASE_URL = "https://www.ncbi.nlm.nih.gov/pmc/articles"

# 搜索策略（多维度覆盖，确保能采集到300篇）
SEARCH_STRATEGIES = [
    # 核心穴位 + 基因表达
    {"query": '("acupuncture"[Title/Abstract] OR "electroacupuncture"[Title/Abstract]) AND ("gene"[Title/Abstract] OR "genomic"[Title/Abstract] OR "transcriptome"[Title/Abstract])', "name": "acupuncture_gene"},
    {"query": '("acupuncture"[Title/Abstract] OR "moxibustion"[Title/Abstract]) AND ("microarray"[Title/Abstract] OR "RNA-seq"[Title/Abstract])', "name": "acupuncture_omics"},
    {"query": '("acupoint"[Title/Abstract] OR "meridian"[Title/Abstract]) AND ("proteomics"[Title/Abstract] OR "metabolomics"[Title/Abstract])', "name": "acupoint_omics"},
    {"query": '"acupuncture"[Title/Abstract] AND ("differential expression"[Title/Abstract] OR "signaling pathway"[Title/Abstract])', "name": "acupuncture_pathway"},
    # 重要穴位名称（英文）
    {"query": '("Zusanli"[Title/Abstract] OR "ST36"[Title/Abstract]) AND ("gene"[Title/Abstract] OR "expression"[Title/Abstract])', "name": "zusanli_gene"},
    {"query": '("Taichong"[Title/Abstract] OR "LR3"[Title/Abstract]) AND ("gene"[Title/Abstract] OR "expression"[Title/Abstract])', "name": "taichong_gene"},
    {"query": '("Baihui"[Title/Abstract] OR "GV20"[Title/Abstract]) AND ("gene"[Title/Abstract] OR "expression"[Title/Abstract])', "name": "baihui_gene"},
    {"query": '("Hegu"[Title/Abstract] OR "LI4"[Title/Abstract]) AND ("gene"[Title/Abstract] OR "expression"[Title/Abstract])', "name": "hegu_gene"},
    {"query": '("Sanyinjiao"[Title/Abstract] OR "SP6"[Title/Abstract]) AND ("gene"[Title/Abstract] OR "expression"[Title/Abstract])', "name": "sanyinjiao_gene"},
    {"query": '("Neiguan"[Title/Abstract] OR "PC6"[Title/Abstract]) AND ("gene"[Title/Abstract] OR "expression"[Title/Abstract])', "name": "neiguan_gene"},
    {"query": '("Fengchi"[Title/Abstract] OR "GB20"[Title/Abstract]) AND ("gene"[Title/Abstract] OR "expression"[Title/Abstract])', "name": "fengchi_gene"},
    {"query": '("Quchi"[Title/Abstract] OR "LI11"[Title/Abstract]) AND ("gene"[Title/Abstract] OR "expression"[Title/Abstract])', "name": "quchi_gene"},
    # 疾病模型 + 穴位
    {"query": '("acupuncture"[Title/Abstract]) AND ("depression"[Title/Abstract] OR "anxiety"[Title/Abstract]) AND ("gene"[Title/Abstract] OR "protein"[Title/Abstract])', "name": "acu_depression_gene"},
    {"query": '("acupuncture"[Title/Abstract]) AND ("hypertension"[Title/Abstract] OR "blood pressure"[Title/Abstract]) AND ("gene"[Title/Abstract] OR "expression"[Title/Abstract])', "name": "acu_hypertension_gene"},
    {"query": '("acupuncture"[Title/Abstract]) AND ("Parkinson"[Title/Abstract] OR "Alzheimer"[Title/Abstract]) AND ("gene"[Title/Abstract] OR "expression"[Title/Abstract])', "name": "acu_neurodegenerative_gene"},
    {"query": '("acupuncture"[Title/Abstract]) AND ("pain"[Title/Abstract] OR "analgesia"[Title/Abstract]) AND ("gene"[Title/Abstract] OR "expression"[Title/Abstract])', "name": "acu_pain_gene"},
    {"query": '("acupuncture"[Title/Abstract]) AND ("inflammation"[Title/Abstract] OR "immune"[Title/Abstract]) AND ("gene"[Title/Abstract] OR "cytokine"[Title/Abstract])', "name": "acu_inflammation_gene"},
    {"query": '("acupuncture"[Title/Abstract]) AND ("cancer"[Title/Abstract] OR "tumor"[Title/Abstract]) AND ("gene"[Title/Abstract] OR "expression"[Title/Abstract])', "name": "acu_cancer_gene"},
    {"query": '("acupuncture"[Title/Abstract]) AND ("diabetes"[Title/Abstract] OR "metabolic"[Title/Abstract]) AND ("gene"[Title/Abstract] OR "expression"[Title/Abstract])', "name": "acu_metabolic_gene"},
    # 机制研究
    {"query": '"acupuncture"[Title/Abstract] AND ("mechanism"[Title/Abstract] OR "molecular"[Title/Abstract]) AND ("network pharmacology"[Title/Abstract] OR "systems biology"[Title/Abstract])', "name": "acu_mechanism_network"},
    {"query": '"acupuncture"[Title/Abstract] AND ("hippocampus"[Title/Abstract] OR "cortex"[Title/Abstract] OR "hypothalamus"[Title/Abstract]) AND ("gene"[Title/Abstract] OR "expression"[Title/Abstract])', "name": "acu_brain_gene"},
    {"query": '"acupuncture"[Title/Abstract] AND ("autophagy"[Title/Abstract] OR "apoptosis"[Title/Abstract] OR "oxidative stress"[Title/Abstract]) AND ("gene"[Title/Abstract] OR "protein"[Title/Abstract])', "name": "acu_cellular_gene"},
]

# 已知的重要穴位列表（用于自动识别）
KNOWN_ACUPOINTS = {
    "英文": [
        "Zusanli", "ST36", "Taichong", "LR3", "Baihui", "GV20", "Hegu", "LI4",
        "Sanyinjiao", "SP6", "Neiguan", "PC6", "Fengchi", "GB20", "Quchi", "LI11",
        "Yanglingquan", "GB34", "Guanyuan", "CV4", "Qihai", "CV6", "Shenshu", "BL23",
        "Pishu", "BL20", "Ganshu", "BL18", "Xinshu", "BL15", "Feishu", "BL13",
        "Dazhui", "GV14", "Mingmen", "GV4", "Changqiang", "GV1", "Yintang", "EX-HN3",
        "Tianshu", "ST25", "Shuigou", "GV26", "Renzhong", "GV26", "Jingming", "BL1",
        "Chengjiang", "CV24", "Jianjing", "GB21", "Fenglong", "ST40", "Kongzui", "LU6",
        "Lieque", "LU7", "Taiyuan", "LU9", "Yujie", "LU10", "Shaoshang", "LU11",
        "Shangyang", "LI1", "Sanjian", "LI3", "Yangxi", "LI5", "Pianli", "LI6",
        "Wenliu", "LI7", "Xialian", "LI8", "Sanli", "LI10", "Jianyu", "LI15",
        "Jianliao", "TE14", "Naohui", "TE13", "Sizhukong", "TE23", "Yifeng", "TE17",
        "Jiaosun", "TE20", "Ermen", "TE21", "Tinggong", "SI19", "Tinghui", "GB2",
        "Shangguan", "GB3", "Xiaguan", "ST7", "Jiache", "ST6", "Dicang", "ST4",
        "Yingxiang", "LI20", "Sibai", "ST2", "Cuanzhu", "BL2", "Meichong", "BL3",
        "Qucha", "BL4", "Wuchu", "BL5", "Chengguang", "BL6", "Tongtian", "BL7",
        "Luoque", "BL8", "Yuzhen", "BL9", "Tianzhu", "BL10", "Dazhu", "BL11",
        "Fengmen", "BL12", "Jueyinshu", "BL14", "Dushu", "BL16", "Geshu", "BL17",
        "Danshu", "BL19", "Weishu", "BL21", "Sanjiaoshu", "BL22", "Dachangshu", "BL25",
        "Guanyuanshu", "BL26", "Xiaochangshu", "BL27", "Pangguangshu", "BL28", "Baihuanshu", "BL30",
        "Shangliao", "BL31", "Ciliao", "BL32", "Zhongliao", "BL33", "Xialiao", "BL34",
        "Huiyang", "BL35", "Fufen", "BL41", "Pohu", "BL42", "Gaohuang", "BL43",
        "Shentang", "BL44", "Yixi", "BL45", "Geguan", "BL46", "Hunmen", "BL47",
        "Yanggang", "BL48", "Yishe", "BL49", "Weicang", "BL50", "Huangmen", "BL51",
        "Zhishi", "BL52", "Baohuang", "BL53", "Zhibian", "BL54", "Heyang", "BL55",
        "Chengjin", "BL56", "Chengshan", "BL57", "Feiyang", "BL58", "Fuyang", "BL59",
        "Kunlun", "BL60", "Pucan", "BL61", "Shenmai", "BL62", "Jinmen", "BL63",
        "Jinggu", "BL64", "Shugu", "BL65", "Zutonggu", "BL66", "Zhiyin", "BL67",
        "Zhiyin", "BL67",
    ],
    "中文": [
        "足三里", "太冲", "百会", "合谷", "三阴交", "内关", "风池", "曲池",
        "阳陵泉", "关元", "气海", "肾俞", "脾俞", "肝俞", "心俞", "肺俞",
        "大椎", "命门", "长强", "印堂", "天枢", "水沟", "人中", "睛明",
        "承浆", "肩井", "丰隆", "孔最", "列缺", "太渊", "鱼际", "少商",
    ]
}

# 请求限速（NCBI无API key: 3 req/sec; 有API key: 10 req/sec）
REQUEST_DELAY = 0.35 if NCBI_API_KEY else 0.35

# 批次大小
BATCH_SIZE = 100
