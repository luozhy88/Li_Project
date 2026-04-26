#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
穴位名称标准化模块
依据《腧穴名称与定位》(GB/T 12346-2021) 国家标准
将各种别名、拼音、英文代码统一为标准中文名称
"""

import sqlite3
import csv
import os
from collections import Counter

# 穴位标准化映射表
# key: 文献中可能出现的名称（大小写不敏感）
# value: 标准中文名称
ACUPOINT_NORMALIZATION_MAP = {
    # 足三里 (ST36)
    'st36': '足三里',
    'zusanli': '足三里',
    'sanli': '足三里',
    
    # 太冲 (LR3)
    'lr3': '太冲',
    'taichong': '太冲',
    
    # 百会 (GV20)
    'gv20': '百会',
    'baihui': '百会',
    
    # 三阴交 (SP6)
    'sp6': '三阴交',
    'sanyinjiao': '三阴交',
    
    # 丰隆 (ST40)
    'st40': '丰隆',
    'fenglong': '丰隆',
    
    # 阳陵泉 (GB34)
    'gb34': '阳陵泉',
    'yanglingquan': '阳陵泉',
    
    # 合谷 (LI4)
    'li4': '合谷',
    'hegu': '合谷',
    
    # 内关 (PC6)
    'pc6': '内关',
    'neiguan': '内关',
    
    # 肾俞 (BL23)
    'bl23': '肾俞',
    'shenshu': '肾俞',
    
    # 大椎 (GV14)
    'gv14': '大椎',
    'dazhui': '大椎',
    
    # 睛明 (BL1)
    'bl1': '睛明',
    'jingming': '睛明',
    
    # 攒竹 (BL2)
    'bl2': '攒竹',
    'cuanzhu': '攒竹',
    
    # 曲差 (BL4)
    'bl4': '曲差',
    'qucha': '曲差',
    
    # 胃俞 (BL21)
    'bl21': '胃俞',
    'weishu': '胃俞',
    
    # 大肠俞 (BL25)
    'bl25': '大肠俞',
    'dachangshu': '大肠俞',
    
    # 天枢 (ST25)
    'st25': '天枢',
    'tianshu': '天枢',
    
    # 地仓 (ST4)
    'st4': '地仓',
    
    # 关元 (CV4 / RN4)
    'cv4': '关元',
    'rn4': '关元',
    'guanyuan': '关元',
    
    # 气海 (CV6 / RN6)
    'cv6': '气海',
    'rn6': '气海',
    'qihai': '气海',
    
    # 脾俞 (BL20)
    'bl20': '脾俞',
    'pishu': '脾俞',
    
    # 肺俞 (BL13)
    'bl13': '肺俞',
    'feishu': '肺俞',
    
    # 肝俞 (BL18)
    'bl18': '肝俞',
    'ganshu': '肝俞',
    
    # 风池 (GB20)
    'gb20': '风池',
    'fengchi': '风池',
    
    # 曲池 (LI11)
    'li11': '曲池',
    'quchi': '曲池',
    
    # 昆仑 (BL60)
    'bl60': '昆仑',
    'kunlun': '昆仑',
    
    # 水沟 (GV26 / DU26)
    'gv26': '水沟',
    'du26': '水沟',
    'shuigou': '水沟',
    'renzhong': '水沟',
    
    # 长强 (GV1)
    'gv1': '长强',
    'changqiang': '长强',
    
    # 命门 (GV4)
    'gv4': '命门',
    'mingmen': '命门',
    
    # 肩髃 (LI15)
    'li15': '肩髃',
    'jianyu': '肩髃',
    
    # 商阳 (LI1)
    'li1': '商阳',
    'shangyang': '商阳',
    
    # 次髎 (BL32)
    'bl32': '次髎',
    'ciliao': '次髎',
    
    # 秩边 (BL54)
    'bl54': '秩边',
    'zhibian': '秩边',
    
    # 印堂 (EX-HN3)
    'yintang': '印堂',
    'ex-hn3': '印堂',
    
    # 眉冲 (BL3)
    'bl3': '眉冲',
    'meichong': '眉冲',
    
    # 五处 (BL5)
    'bl5': '五处',
    'wuchu': '五处',
    
    # 承光 (BL6)
    'bl6': '承光',
    'chengguang': '承光',
    
    # 听会 (GB2)
    'gb2': '听会',
    'tinghui': '听会',
    
    # 上关 (GB3)
    'gb3': '上关',
    'shangguan': '上关',
    
    # 四白 (ST2)
    'st2': '四白',
    'sibai': '四白',
    
    # 耳门 (TE21 / SJ21)
    'ermen': '耳门',
    'te21': '耳门',
    'sj21': '耳门',
    
    # 大杼 (BL11) - 注意Dazhu可能是大杼或大椎的误写
    'bl11': '大杼',
    'dazhu': '大杼',
    
    # 其他可能的标准化名称（已有中文的直接保留）
    '足三里': '足三里',
    '太冲': '太冲',
    '百会': '百会',
    '三阴交': '三阴交',
    '丰隆': '丰隆',
    '阳陵泉': '阳陵泉',
    '合谷': '合谷',
    '内关': '内关',
    '肾俞': '肾俞',
    '大椎': '大椎',
    '睛明': '睛明',
    '天枢': '天枢',
    '地仓': '地仓',
    '关元': '关元',
    '气海': '气海',
    '脾俞': '脾俞',
    '肺俞': '肺俞',
    '肝俞': '肝俞',
    '风池': '风池',
    '曲池': '曲池',
    '昆仑': '昆仑',
    '水沟': '水沟',
    '长强': '长强',
    '命门': '命门',
    '肩髃': '肩髃',
    '商阳': '商阳',
    '次髎': '次髎',
    '秩边': '秩边',
    '印堂': '印堂',
    '眉冲': '眉冲',
    '五处': '五处',
    '承光': '承光',
    '听会': '听会',
    '上关': '上关',
    '四白': '四白',
    '耳门': '耳门',
    '大杼': '大杼',
}

# 经络归属映射
MERIDIAN_MAP = {
    '足三里': '足阳明胃经',
    '太冲': '足厥阴肝经',
    '百会': '督脉',
    '三阴交': '足太阴脾经',
    '丰隆': '足阳明胃经',
    '阳陵泉': '足少阳胆经',
    '合谷': '手阳明大肠经',
    '内关': '手厥阴心包经',
    '肾俞': '足太阳膀胱经',
    '大椎': '督脉',
    '睛明': '足太阳膀胱经',
    '攒竹': '足太阳膀胱经',
    '曲差': '足太阳膀胱经',
    '胃俞': '足太阳膀胱经',
    '大肠俞': '足太阳膀胱经',
    '天枢': '足阳明胃经',
    '地仓': '足阳明胃经',
    '关元': '任脉',
    '气海': '任脉',
    '脾俞': '足太阳膀胱经',
    '肺俞': '足太阳膀胱经',
    '肝俞': '足太阳膀胱经',
    '风池': '足少阳胆经',
    '曲池': '手阳明大肠经',
    '昆仑': '足太阳膀胱经',
    '水沟': '督脉',
    '长强': '督脉',
    '命门': '督脉',
    '肩髃': '手阳明大肠经',
    '商阳': '手阳明大肠经',
    '次髎': '足太阳膀胱经',
    '秩边': '足太阳膀胱经',
    '印堂': '经外奇穴',
    '眉冲': '足太阳膀胱经',
    '五处': '足太阳膀胱经',
    '承光': '足太阳膀胱经',
    '听会': '足少阳胆经',
    '上关': '足少阳胆经',
    '四白': '足阳明胃经',
    '耳门': '手少阳三焦经',
    '大杼': '足太阳膀胱经',
    '攒竹': '足太阳膀胱经',
    '曲差': '足太阳膀胱经',
    '胃俞': '足太阳膀胱经',
    '大肠俞': '足太阳膀胱经',
}

# 穴位代码映射
CODE_MAP = {
    '足三里': 'ST36',
    '太冲': 'LR3',
    '百会': 'GV20',
    '三阴交': 'SP6',
    '丰隆': 'ST40',
    '阳陵泉': 'GB34',
    '合谷': 'LI4',
    '内关': 'PC6',
    '肾俞': 'BL23',
    '大椎': 'GV14',
    '睛明': 'BL1',
    '攒竹': 'BL2',
    '曲差': 'BL4',
    '胃俞': 'BL21',
    '大肠俞': 'BL25',
    '天枢': 'ST25',
    '地仓': 'ST4',
    '关元': 'CV4',
    '气海': 'CV6',
    '脾俞': 'BL20',
    '肺俞': 'BL13',
    '肝俞': 'BL18',
    '风池': 'GB20',
    '曲池': 'LI11',
    '昆仑': 'BL60',
    '水沟': 'GV26',
    '长强': 'GV1',
    '命门': 'GV4',
    '肩髃': 'LI15',
    '商阳': 'LI1',
    '次髎': 'BL32',
    '秩边': 'BL54',
    '印堂': 'EX-HN3',
    '眉冲': 'BL3',
    '五处': 'BL5',
    '承光': 'BL6',
    '听会': 'GB2',
    '上关': 'GB3',
    '四白': 'ST2',
    '耳门': 'TE21',
    '大杼': 'BL11',
    '攒竹': 'BL2',
    '曲差': 'BL4',
    '胃俞': 'BL21',
    '大肠俞': 'BL25',
}


def normalize_acupoint_name(name):
    """将穴位名称标准化为标准中文名称"""
    if not name:
        return None
    name_lower = name.strip().lower()
    normalized = ACUPOINT_NORMALIZATION_MAP.get(name_lower)
    if normalized:
        return normalized
    # 如果已经是中文标准名，直接返回
    if name in MERIDIAN_MAP:
        return name
    # 未匹配到的，打印警告并返回原名
    print(f"[警告] 未找到穴位标准化映射: '{name}'，保留原名")
    return name


def normalize_database(db_path='output/acupoint_gene.db'):
    """标准化数据库中的穴位名称"""
    if not os.path.exists(db_path):
        print(f"数据库文件不存在: {db_path}")
        return
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # 1. 获取所有穴位
    cursor.execute('SELECT id, name_en, name_cn FROM acupoints')
    acupoints = cursor.fetchall()
    
    # 建立标准化映射
    name_mapping = {}  # old_name -> normalized_name
    id_mapping = {}    # id -> normalized_name
    
    for aid, name_en, name_cn in acupoints:
        old_name = name_en or name_cn or ''
        normalized = normalize_acupoint_name(old_name)
        if normalized:
            name_mapping[old_name] = normalized
            id_mapping[aid] = normalized
    
    # 2. 检查是否有冲突（多个不同旧名映射到同一个标准化名称）
    normalized_to_ids = {}
    for aid, norm_name in id_mapping.items():
        normalized_to_ids.setdefault(norm_name, []).append(aid)
    
    # 3. 创建新的标准化acupoints表
    # 先删除旧表数据，插入标准化后的数据
    cursor.execute('DELETE FROM acupoints')
    
    seen_names = set()
    new_id_map = {}  # normalized_name -> new_id
    
    for norm_name, ids in normalized_to_ids.items():
        if norm_name in seen_names:
            continue
        seen_names.add(norm_name)
        
        code = CODE_MAP.get(norm_name, '')
        meridian = MERIDIAN_MAP.get(norm_name, '')
        
        cursor.execute(
            'INSERT INTO acupoints (name_en, name_cn, code, meridian) VALUES (?, ?, ?, ?)',
            (norm_name, norm_name, code, meridian)
        )
        new_id = cursor.lastrowid
        new_id_map[norm_name] = new_id
    
    # 4. 更新 article_acupoint 关联表
    cursor.execute('SELECT article_id, acupoint_id FROM article_acupoint')
    article_acupoint_rows = cursor.fetchall()
    
    cursor.execute('DELETE FROM article_acupoint')
    
    # 统计每对 article-acupoint 出现次数（去重）
    unique_pairs = set()
    for article_id, acupoint_id in article_acupoint_rows:
        norm_name = id_mapping.get(acupoint_id)
        if norm_name and norm_name in new_id_map:
            new_acupoint_id = new_id_map[norm_name]
            unique_pairs.add((article_id, new_acupoint_id))
    
    for article_id, acupoint_id in unique_pairs:
        cursor.execute(
            'INSERT INTO article_acupoint (article_id, acupoint_id) VALUES (?, ?)',
            (article_id, acupoint_id)
        )
    
    conn.commit()
    conn.close()
    
    print(f"数据库穴位标准化完成！")
    print(f"标准化后穴位数量: {len(seen_names)}")
    for name in sorted(seen_names):
        print(f"  - {name} ({CODE_MAP.get(name, '?')})")


def update_acupoints_csv(csv_path='output/acupoints.csv'):
    """更新CSV文件中的穴位名称"""
    if not os.path.exists(csv_path):
        print(f"CSV文件不存在: {csv_path}")
        return
    
    # 读取原CSV，统计各穴位出现次数
    counter = Counter()
    with open(csv_path, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row.get('name_en', '') or row.get('name_cn', '')
            normalized = normalize_acupoint_name(name)
            if normalized:
                counter[normalized] += int(row.get('article_count', 1))
    
    # 写回标准化后的CSV
    with open(csv_path, 'w', encoding='utf-8-sig', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['name_cn', 'code', 'meridian', 'article_count'])
        for name, count in counter.most_common():
            writer.writerow([name, CODE_MAP.get(name, ''), MERIDIAN_MAP.get(name, ''), count])
    
    print(f"CSV文件已更新: {csv_path}")


def generate_summary():
    """生成标准化后的统计摘要"""
    conn = sqlite3.connect('output/acupoint_gene.db')
    cursor = conn.cursor()
    
    cursor.execute('''
        SELECT a.name_cn, a.code, a.meridian, COUNT(DISTINCT aa.article_id) as article_count
        FROM acupoints a
        LEFT JOIN article_acupoint aa ON a.id = aa.acupoint_id
        GROUP BY a.id
        ORDER BY article_count DESC
    ''')
    
    results = cursor.fetchall()
    conn.close()
    
    print("\n=== 标准化后穴位统计 ===")
    print(f"{'穴位名称':<10} {'代码':<8} {'经络':<16} {'文献数':<8}")
    print("-" * 50)
    for name, code, meridian, count in results:
        print(f"{name:<10} {code:<8} {meridian:<16} {count:<8}")
    
    return results


if __name__ == '__main__':
    normalize_database()
    update_acupoints_csv()
    generate_summary()
