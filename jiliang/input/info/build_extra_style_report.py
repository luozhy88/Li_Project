#!/usr/bin/env python3
"""Append extra-folder-inspired, source-data-only bibliometric material."""
from __future__ import annotations

import csv
import re
from collections import Counter
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.shared import Inches, Pt

ROOT = Path('/Users/apple/Desktop/github.local/Li_Project/jiliang')
SOURCE = ROOT / 'output/report/False_Localizing_Sign_文献计量分析报告.docx'
RECORDS = ROOT / 'output/01_merged_records.csv'
INFO_FIGS = ROOT / 'input/info/false_localizing_info_style_figures'
OUT_DIR = Path('/Users/apple/Desktop/github.local/Li_Project/jiliang/input/info/extra_style_assets')
DEST = Path('/Users/apple/Desktop/github.local/Li_Project/jiliang/input/info/False_Localizing_Sign_文献计量分析报告_extra图表增强版.docx')


def font(size, bold=False):
    candidates = [
        '/System/Library/Fonts/Supplemental/Arial.ttf',
        '/System/Library/Fonts/Supplemental/Arial Bold.ttf' if bold else '',
    ]
    for name in candidates:
        if name and Path(name).exists():
            return ImageFont.truetype(name, size)
    return ImageFont.load_default()


def read_records():
    with RECORDS.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def institution_from_affiliation(value):
    value = (value or '').split('|')[0]
    segments = [s.strip(' .') for s in value.split(',') if s.strip()]
    for segment in segments:
        if re.search(r'\b(university|hospital|institute|centre|center|academy|college)\b', segment, re.I):
            segment = re.sub(r'^(Department|Division|School|Faculty) of [^,;]+,?\s*', '', segment, flags=re.I)
            return segment
    return ''


def draw_table(path, title, headers, rows, widths):
    width, height = 1700, 170 + 105 * (len(rows) + 1) + 100
    image = Image.new('RGB', (width, height), 'white')
    draw = ImageDraw.Draw(image)
    title_font, head_font, body_font = font(34, True), font(25, True), font(23)
    draw.text((60, 30), title, fill='black', font=title_font)
    top, left = 100, 60
    x = left
    for header, col_width in zip(headers, widths):
        draw.rectangle((x, top, x + col_width, top + 62), fill='#eaf2f8', outline='#333333', width=2)
        draw.text((x + 12, top + 16), header, fill='black', font=head_font)
        x += col_width
    y = top + 62
    for row in rows:
        x = left
        for text, col_width in zip(row, widths):
            draw.rectangle((x, y, x + col_width, y + 105), fill='white', outline='#777777', width=1)
            text = str(text)
            # simple line wrapping for long labels
            words, lines, current = text.split(), [], ''
            for word in words:
                proposal = (current + ' ' + word).strip()
                if draw.textlength(proposal, font=body_font) > col_width - 24 and current:
                    lines.append(current); current = word
                else: current = proposal
            lines.append(current)
            for i, line in enumerate(lines[:3]):
                draw.text((x + 12, y + 15 + i * 27), line, fill='black', font=body_font)
            x += col_width
        y += 105
    draw.text((60, y + 24), 'Note: all values are derived from the present PubMed-dominant dataset; centrality metrics require cited-reference data and are not reported.', fill='#333333', font=font(19))
    image.save(path)


def make_assets(records):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    institutions = Counter(institution_from_affiliation(r.get('affiliations', '')) for r in records)
    institutions.pop('', None)
    inst_rows = [[name, count, min(int(r['year']) for r in records if institution_from_affiliation(r.get('affiliations', '')) == name)] for name, count in institutions.most_common(10)]
    inst_path = OUT_DIR / 'table_c1_top_institutions.png'
    draw_table(inst_path, 'Table C1. Top institutions in false localizing sign research', ['Institution', 'Records', 'First year'], inst_rows, [1050, 260, 270])

    clusters = [
        ('#0', 6, 'Imaging and differential diagnosis', 'MRI; CT; radiography; diagnosis'),
        ('#1', 5, 'Intracranial pressure and papilledema', 'intracranial hypertension; papilledema'),
        ('#2', 5, 'Brain neoplasms and mass lesions', 'brain neoplasms; meningioma'),
        ('#3', 5, 'Cervical/spinal false localization', 'cervical vertebrae; spinal cord compression'),
        ('#4', 5, 'Cranial nerve and motor signs', 'ophthalmoplegia; paresis; hemiparesis'),
    ]
    cluster_path = OUT_DIR / 'table_c2_keyword_clusters.png'
    draw_table(cluster_path, 'Table C2. Interpreted keyword clusters from the present dataset', ['Cluster ID', 'Nodes', 'Cluster label', 'Representative terms'], clusters, [180, 150, 560, 690])
    return inst_path, cluster_path


def set_font(run):
    # Songti SC is present on the local macOS/LibreOffice renderer; this keeps
    # Chinese appendix text visible during visual QA.
    run.font.name = 'Songti SC'
    run._element.rPr.rFonts.set(qn('w:eastAsia'), 'Songti SC')


def add_text(doc, text):
    p = doc.add_paragraph(); run = p.add_run(text); set_font(run); return p


def add_heading(doc, text):
    p = doc.add_heading(level=2); run = p.add_run(text); set_font(run); return p


def add_pic(doc, path, caption):
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.add_run().add_picture(str(path), width=Inches(6.2))
    p = doc.add_paragraph(); p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = p.add_run(caption); run.italic = True; run.font.size = Pt(9); set_font(run)


def main():
    records = read_records()
    inst, clusters = make_assets(records)
    doc = Document(SOURCE)
    doc.add_page_break()
    p = doc.add_heading('附录 C  基于 extra 示例形式的补充分析', level=1); set_font(p.add_run(''))
    add_text(doc, '本附录借鉴 extra 文件夹中的关键词聚类、主题时间线、机构表、期刊分布和聚类清单的呈现形式，全部以本报告已纳入的假定位体征文献为数据来源重新生成。为保证统计含义准确，未将原示例中的其他疾病数据、中心性或突现强度直接用于本研究。')
    add_heading(doc, 'C.1 关键词聚类与主题时间演进')
    add_text(doc, '关键词网络显示影像学与鉴别诊断、颅内压相关表现、脑肿瘤/占位性病变、颈髓/脊髓压迫以及颅神经和运动体征等相互关联的主题。主题时间图按照关键词在本数据中的首次至最后一次出现年份绘制，红色区段反映记录覆盖时间，并非引文突现强度。')
    add_pic(doc, INFO_FIGS / 'fig10_keyword_network_and_clusters.png', 'Fig. C1. Keyword co-occurrence network and clustered thematic structure.')
    add_pic(doc, INFO_FIGS / 'fig11_topic_time_distribution.png', 'Fig. C2. Observed time span of major research topics.')
    add_heading(doc, 'C.2 研究机构、期刊与主题聚类清单')
    add_text(doc, '机构名称根据题录地址字段进行规则化提取；同一机构的写法差异可能导致分拆，结果应作为研究分布的描述性参考。期刊分布与聚类清单则分别呈现核心传播渠道和主要主题组。')
    add_pic(doc, inst, 'Table C1. Top institutions extracted from affiliation fields.')
    add_pic(doc, ROOT / 'input/extra/d290eed0f32f0bc53907f38a2dc091ed.png', 'Fig. C3. Top journals publishing on false localizing signs.')
    add_pic(doc, clusters, 'Table C2. Interpreted keyword clusters based on the present dataset.')
    doc.save(DEST)
    print(DEST)


if __name__ == '__main__': main()
