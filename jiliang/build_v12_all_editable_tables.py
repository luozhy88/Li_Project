#!/usr/bin/env python3
"""Create v12 with all report tables as editable three-line Word tables."""
from __future__ import annotations

import csv
import json
import re
from collections import Counter
from pathlib import Path

from docx import Document
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt


ROOT = Path('/Users/apple/Desktop/github.local/Li_Project/jiliang')
SOURCE = ROOT / 'output/report/10_False_Localizing_Sign_文献计量分析报告_Medicine投稿修订版v11.docx'
DEST = ROOT / 'output/report/10_False_Localizing_Sign_文献计量分析报告_Medicine投稿修订版v12.docx'
DATA = ROOT / 'output/01_merged_records.csv'
STATS = ROOT / 'output/logs/b1_own_stats.json'


def set_cell_margins(cell, top=80, start=120, bottom=80, end=120):
    tc_pr = cell._tc.get_or_add_tcPr()
    node = tc_pr.first_child_found_in('w:tcMar')
    if node is None:
        node = OxmlElement('w:tcMar'); tc_pr.append(node)
    for side, value in [('top', top), ('start', start), ('bottom', bottom), ('end', end)]:
        side_node = node.find(qn(f'w:{side}'))
        if side_node is None:
            side_node = OxmlElement(f'w:{side}'); node.append(side_node)
        side_node.set(qn('w:w'), str(value)); side_node.set(qn('w:type'), 'dxa')


def set_cell_width(cell, width):
    tc_pr = cell._tc.get_or_add_tcPr()
    node = tc_pr.first_child_found_in('w:tcW')
    if node is None:
        node = OxmlElement('w:tcW'); tc_pr.append(node)
    node.set(qn('w:w'), str(width)); node.set(qn('w:type'), 'dxa')


def set_geometry(table, widths):
    table.autofit = False
    pr = table._tbl.tblPr
    tbl_w = pr.first_child_found_in('w:tblW')
    if tbl_w is None:
        tbl_w = OxmlElement('w:tblW'); pr.insert(0, tbl_w)
    tbl_w.set(qn('w:w'), str(sum(widths))); tbl_w.set(qn('w:type'), 'dxa')
    ind = pr.first_child_found_in('w:tblInd')
    if ind is None:
        ind = OxmlElement('w:tblInd'); pr.append(ind)
    ind.set(qn('w:w'), '120'); ind.set(qn('w:type'), 'dxa')
    layout = pr.first_child_found_in('w:tblLayout')
    if layout is None:
        layout = OxmlElement('w:tblLayout'); pr.append(layout)
    layout.set(qn('w:type'), 'fixed')
    grid = table._tbl.tblGrid
    for child in list(grid): grid.remove(child)
    for width in widths:
        col = OxmlElement('w:gridCol'); col.set(qn('w:w'), str(width)); grid.append(col)
    for row in table.rows:
        for cell, width in zip(row.cells, widths):
            set_cell_width(cell, width); set_cell_margins(cell)


def set_three_lines(table):
    pr = table._tbl.tblPr
    borders = pr.first_child_found_in('w:tblBorders')
    if borders is None:
        borders = OxmlElement('w:tblBorders'); pr.append(borders)
    for edge in ['top', 'left', 'bottom', 'right', 'insideH', 'insideV']:
        node = borders.find(qn(f'w:{edge}'))
        if node is None:
            node = OxmlElement(f'w:{edge}'); borders.append(node)
        node.set(qn('w:val'), 'nil')
    for edge, size in [('top', '12'), ('bottom', '12')]:
        node = borders.find(qn(f'w:{edge}'))
        node.set(qn('w:val'), 'single'); node.set(qn('w:sz'), size); node.set(qn('w:color'), '000000')
    for cell in table.rows[0].cells:
        tc_pr = cell._tc.get_or_add_tcPr()
        cb = tc_pr.first_child_found_in('w:tcBorders')
        if cb is None:
            cb = OxmlElement('w:tcBorders'); tc_pr.append(cb)
        node = OxmlElement('w:bottom')
        node.set(qn('w:val'), 'single'); node.set(qn('w:sz'), '8'); node.set(qn('w:color'), '000000')
        cb.append(node)


def put_text(cell, text, bold=False, center=False, size=9):
    cell.text = ''; cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
    p = cell.paragraphs[0]
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER if center else WD_ALIGN_PARAGRAPH.LEFT
    p.paragraph_format.space_before = Pt(0); p.paragraph_format.space_after = Pt(0); p.paragraph_format.line_spacing = 1.05
    run = p.add_run(str(text)); run.bold = bold; run.font.name = 'Times New Roman'; run.font.size = Pt(size)
    run._element.rPr.rFonts.set(qn('w:ascii'), 'Times New Roman'); run._element.rPr.rFonts.set(qn('w:hAnsi'), 'Times New Roman')


def add_table(doc, headers, rows, widths, sizes=None):
    table = doc.add_table(rows=1, cols=len(headers)); table.alignment = WD_TABLE_ALIGNMENT.LEFT
    set_geometry(table, widths)
    for cell, header in zip(table.rows[0].cells, headers): put_text(cell, header, bold=True, center=True, size=9)
    for row in rows:
        cells = table.add_row().cells
        for i, (cell, value) in enumerate(zip(cells, row)):
            center = i < 2 and len(headers) >= 3
            put_text(cell, value, center=center, size=(sizes[i] if sizes else 9))
    set_geometry(table, widths); set_three_lines(table)
    return table


def replace_image_with_table(doc, caption_prefix, headers, rows, widths, sizes=None, note=None):
    caption = next(p for p in doc.paragraphs if p.text.strip().startswith(caption_prefix))
    old = caption._p.getprevious(); drawing = qn('w:drawing')
    while old is not None and not old.findall('.//' + drawing): old = old.getprevious()
    if old is None: raise RuntimeError(f'Image for {caption_prefix} not found')
    parent = caption._p.getparent(); parent.remove(old)
    caption.text = caption_prefix
    for run in caption.runs:
        run.font.name = 'Times New Roman'; run.font.size = Pt(10); run.italic = True
    table = add_table(doc, headers, rows, widths, sizes)
    tbl = table._tbl; parent.remove(tbl)
    parent.insert(parent.index(caption._p), caption._p)
    parent.insert(parent.index(caption._p) + 1, tbl)
    if note:
        p = doc.add_paragraph(); p.paragraph_format.space_before = Pt(3); p.paragraph_format.space_after = Pt(6)
        r = p.add_run(note); r.font.name = 'Times New Roman'; r.font.size = Pt(8.5); r.italic = True
        parent.remove(p._p); parent.insert(parent.index(tbl) + 1, p._p)


def institution_from_affiliation(value):
    value = (value or '').split('|')[0]
    segments = [s.strip(' .') for s in value.split(',') if s.strip()]
    for segment in segments:
        if re.search(r'\b(university|hospital|institute|centre|center|academy|college)\b', segment, re.I):
            return re.sub(r'^(Department|Division|School|Faculty) of [^,;]+,?\s*', '', segment, flags=re.I)
    return ''


def main():
    records = list(csv.DictReader(DATA.open(encoding='utf-8-sig', newline='')))
    doc = Document(SOURCE)

    inst = Counter(institution_from_affiliation(r.get('affiliations', '')) for r in records)
    inst.pop('', None)
    inst_rows = []
    for name, count in inst.most_common(10):
        years = [int(r['year']) for r in records if institution_from_affiliation(r.get('affiliations', '')) == name and r.get('year')]
        inst_rows.append([name, count, min(years) if years else ''])
    replace_image_with_table(
        doc, 'Table B1. Top institutions extracted from affiliation fields.',
        ['Institution', 'Records', 'First year'], inst_rows, [5400, 1500, 2460], [8.5, 9, 9],
        'Note: institution names were rule-extracted from affiliation fields; spelling variants may be split.')

    stats = json.loads(STATS.read_text(encoding='utf-8'))
    # Table B2 is already a native table in v11; its content is sourced from the same stats JSON as Fig. B2.

    first = Counter(r.get('first_author', '').strip() for r in records if r.get('first_author', '').strip())
    first_rows = []
    for author, count in first.most_common(10):
        own = [r for r in records if r.get('first_author', '').strip() == author]
        chosen = sorted(own, key=lambda r: int(r.get('year') or 0), reverse=True)[0]
        title = ' '.join((chosen.get('title') or '').split())
        first_rows.append([author, f'{count} ({count / len(records) * 100:.2f}%)', title])
    replace_image_with_table(
        doc, 'Table C1. Top 10 productive first authors in false localizing sign research.',
        ['First author', 'Documents', 'Representative first-authored publication'], first_rows,
        [1700, 1500, 6160], [8.5, 8.5, 8.0],
        f'Note: each of the {len(records)} records is credited to its first author only, so every publication appears exactly once.')

    doc.save(DEST)
    print(DEST)


if __name__ == '__main__': main()
