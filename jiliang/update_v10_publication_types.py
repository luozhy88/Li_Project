from pathlib import Path
import shutil
import zipfile
import xml.etree.ElementTree as ET

from PIL import Image, ImageDraw, ImageFont


ROOT = Path('/Users/apple/Desktop/github.local/Li_Project/jiliang')
SOURCE = ROOT / 'output/report/10_False_Localizing_Sign_文献计量分析报告_Medicine投稿修订版v9.docx'
DEST = ROOT / 'output/report/10_False_Localizing_Sign_文献计量分析报告_Medicine投稿修订版v10.docx'
CHART = ROOT / 'output/figures/geo_style/fig_e4_publication_types_v10.png'


def font(size, bold=False):
    path = '/System/Library/Fonts/Supplemental/Arial Bold.ttf' if bold else '/System/Library/Fonts/Supplemental/Arial.ttf'
    return ImageFont.truetype(path, size)


def draw_chart(path):
    """Draw the six-category publication-type chart at the original pixel size."""
    W, H = 2079, 1100
    im = Image.new('RGB', (W, H), 'white')
    d = ImageDraw.Draw(im)
    title = font(42, True)
    label = font(30)
    small = font(28)
    d.text((W // 2, 38), 'Publication types of the included records', anchor='ma', fill='#222222', font=title)
    left, right = 520, 1830
    top, bottom = 175, 960
    cats = [
        ('Case Reports', 72, '#4d9ed7'),
        ('Review (incl. Systematic Review)', 11, '#e58d47'),
        ('Clinical Trial', 1, '#47a85c'),
        ('Historical / Biography', 4, '#a46ccc'),
        ('Journal Article (other)', 66, '#d2af34'),
        ('Other', 1, '#e85d5d'),
    ]
    maxv = 84
    row_h = (bottom - top) / len(cats)
    bar_h = 52
    for i, (name, value, color) in enumerate(cats):
        y = top + row_h * i + row_h / 2
        d.text((left - 24, y), name, anchor='rm', fill='#222222', font=label)
        x2 = left + (right - left) * value / maxv
        d.rectangle((left, y - bar_h / 2, x2, y + bar_h / 2), fill=color)
        d.text((x2 + 16, y), f'{value}  ({value / 155 * 100:.1f}%)', anchor='lm', fill='#333333', font=small)
    d.text(((left + right) / 2, H - 55), 'Number of publications', anchor='ma', fill='#222222', font=label)
    im.save(path)


def main():
    # Exclude the four Letter / Editorial / Comment records from the figure's
    # comparative publication-type denominator.
    type_counts = [
        ('Case Reports', 72),
        ('Review (incl. Systematic Review)', 11),
        ('Clinical Trial', 1),
        ('Historical / Biography', 4),
        ('Journal Article (other)', 66),
        ('Other', 1),
    ]
    chart_total = 155
    CHART.parent.mkdir(parents=True, exist_ok=True)
    draw_chart(CHART)

    shutil.copy2(SOURCE, DEST)

    body_text = (
        '经检索，PubMed 获得 125 条记录，WoS 获得 89 条记录；依次按 PMID、DOI 与归一化标题去重，'
        '共移除 55 条重复记录，最终纳入 159 篇文献进行分析（PubMed 来源 123 篇，WoS 来源 36 篇），'
        '时间跨度为 1904—2026 年。其中具有 DOI 的记录 135 条，具有摘要的记录 123 条。'
        '文献体裁构成按优先级进行互斥归类。由于来信、述评或评论类记录 4 篇不纳入本图的主体出版类型比较，'
        '图 2 基于其余 155 条记录绘制：病例报告 72 篇（46.5%）、综述（含系统综述）11 篇（7.1%）、'
        '临床试验 1 篇（0.6%）、历史/传记类 4 篇（2.6%）、其他期刊论文 66 篇（42.6%）及其他类型 1 篇（0.6%）。'
    )
    caption_text = (
        'Fig. 2. Distribution of publication types among the 155 included records after excluding '
        'Letter/Editorial/Comment records (exclusive single-label classification).'
    )

    ns = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
    w_t = '{%s}t' % ns['w']
    with zipfile.ZipFile(DEST, 'r') as zin:
        files = {name: zin.read(name) for name in zin.namelist()}
    root = ET.fromstring(files['word/document.xml'])
    def replace_para(p, text):
        ppr = p.find('{%s}pPr' % ns['w'])
        for child in list(p):
            if child is not ppr:
                p.remove(child)
        r = ET.SubElement(p, '{%s}r' % ns['w'])
        ET.SubElement(r, w_t).text = text

    for p in root.findall('.//{%s}body/{%s}p' % (ns['w'], ns['w'])):
        text = ''.join(x.text or '' for x in p.iter(w_t))
        if '具有摘要的记录' in text:
            replace_para(p, body_text)
        elif text.strip().startswith('Fig. 2.'):
            replace_para(p, caption_text)
    files['word/document.xml'] = ET.tostring(root, encoding='utf-8', xml_declaration=True)
    files['word/media/image16.png'] = CHART.read_bytes()
    with zipfile.ZipFile(DEST, 'w', zipfile.ZIP_DEFLATED) as zout:
        for name, data in files.items():
            zout.writestr(name, data)
    print(DEST)


if __name__ == '__main__':
    main()
