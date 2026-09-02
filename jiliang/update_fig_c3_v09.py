#!/usr/bin/env python3
"""Update report v08 -> v09: swap in the reviewer-revised Fig. C3 (country
collaboration network) and update its body text and caption.

Fig. C3 fixes: legend added; decorative halo removed; node AREA strictly
proportional to publication output; hand-placed layout (no label overlap, no
edge crossings); isolated nodes aligned in a bottom strip.
"""
from __future__ import annotations

from pathlib import Path

from docx import Document
from docx.oxml.ns import qn
from docx.shared import Inches
from docx.text.paragraph import Paragraph

ROOT = Path('/Users/apple/Desktop/github.local/Li_Project/jiliang')
SOURCE = ROOT / 'output/report/08_False_Localizing_Sign_文献计量分析报告_Medicine投稿修订版v2.docx'
DEST = ROOT / 'output/report/09_False_Localizing_Sign_文献计量分析报告_Medicine投稿修订版v3.docx'
IMG = ROOT / 'output/figures/geo_style/fig_e2_country_network_glow.png'

BODY = ('图 C3 展示国家/地区间的合作网络：节点面积与发文量严格成正比（见图内 Legend 示例），'
        '节点采用黄—橙—红多层配色，红色连线表示两国/地区在同一篇文献中共同署名；'
        '无合作连线的国家/地区（孤立节点）统一排列于图谱底部区域，'
        '中央画幅保留给具有合作关系的国家/地区聚类。')

CAPTION = ('Fig. C3. Network visualization of international collaboration among countries in '
           'false localizing sign research. Each node represents a country, with node area '
           'proportional to publication output; links indicate co-authorship within one record. '
           'Countries without co-authorship links (isolated nodes) are arranged in the bottom '
           'strip. Node sizes follow the scaling shown in the legend.')


def rewrite(p, text):
    if p.runs:
        p.runs[0].text = text
        for r in p.runs[1:]:
            r._element.getparent().remove(r._element)
    else:
        p.add_run(text)


def main():
    doc = Document(SOURCE)
    for p in doc.paragraphs:
        t = p.text.strip()
        if t.startswith('图 C3 展示国家/地区间的合作网络'):
            rewrite(p, BODY)
            print('正文 C.2 段落已更新')
        elif t.startswith('Fig. C3.'):
            el = p._p.getprevious()
            while el is not None and not el.findall('.//' + qn('w:drawing')):
                el = el.getprevious()
            assert el is not None, '未找到 Fig. C3 对应的图片段落'
            pic_par = Paragraph(el, p._parent)
            pic_par.clear()
            pic_par.alignment = p.alignment  # 保持居中
            pic_par.add_run().add_picture(str(IMG), width=Inches(6.25))
            rewrite(p, CAPTION)
            print('Fig. C3 图片与图注已更新')
    doc.save(DEST)
    print(DEST)


if __name__ == '__main__':
    main()
