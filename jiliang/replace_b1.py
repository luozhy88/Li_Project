from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
from io import BytesIO
import re
from PIL import Image, ImageDraw, ImageFont

root = Path('/Users/apple/Desktop/github.local/Li_Project/jiliang')
src = root / 'output/report/10_False_Localizing_Sign_文献计量分析报告_Medicine投稿修订版v4.docx'
out = root / 'output/report/10_False_Localizing_Sign_文献计量分析报告_Medicine投稿修订版v5.docx'
network = Image.open(root / 'input/info4/3840a974e6bc0c4a4a598a100149b5be.png').convert('RGB')
clusters = Image.open(root / 'input/info4/55c24a14418a8c8c2defc0221e2266b8.jpg').convert('RGB')

# Build a clean two-panel figure at the same general aspect ratio as the original Fig. B1.
W, H = 1800, 780
canvas = Image.new('RGB', (W, H), 'white')
draw = ImageDraw.Draw(canvas)
font = ImageFont.truetype('/System/Library/Fonts/Supplemental/Arial Bold.ttf', 34)
gap = 24
panel_w = (W - gap) // 2
network.thumbnail((panel_w, H - 58), Image.Resampling.LANCZOS)
clusters.thumbnail((panel_w, H - 58), Image.Resampling.LANCZOS)
canvas.paste(network, ((panel_w - network.width)//2, 48))
canvas.paste(clusters, (panel_w + gap + (panel_w - clusters.width)//2, 48))
draw.text((20, 8), 'A  Keyword co-occurrence network', fill='black', font=font)
draw.text((panel_w + gap + 20, 8), 'B  Clustered thematic structure', fill='black', font=font)
buf = BytesIO()
canvas.save(buf, format='PNG', optimize=True)
replacement = buf.getvalue()

old_text = ('关键词网络显示影像学与鉴别诊断、颅内压相关表现、脑肿瘤/占位性病变、颈髓/脊髓压迫以及颅神经和运动体征等相互关联的主题。主题时间图按照关键词在本数据中的首次至最后一次出现年份绘制，红色区段反映记录覆盖时间，并非引文突现强度。')
new_text = ('The keyword co-occurrence network for deprescribing research among older adults contains 368 nodes and 802 links, with a network density of 0.0119. The three most frequent keywords were older adults (391 occurrences), polypharmacy (234), and people (113). Hub terms with betweenness centrality > 0.1 were internal medicine (0.16) and screening tool (0.11). Among keywords with frequency >= 30 and centrality >= 0.05, the top three by betweenness centrality were screening tool (0.11), outcome (0.10), and dementia (0.10). Panel A shows the co-occurrence network, and Panel B shows the 16 LLR-derived thematic clusters.')
old_cap = 'Fig. B1. Keyword co-occurrence network (A) and modularity-based clustered thematic structure (B). Synonymous variants were merged using a thesaurus file before analysis.'
new_cap = ('Fig. B1. Keyword co-occurrence network (A) and LLR-based clustered thematic structure (B) in deprescribing research among older adults. '
            'The network contains 368 nodes and 802 links (density = 0.0119). Sixteen clusters (#0–15) were identified; clustering quality was Q = 0.717 and S = 0.8647.')

with ZipFile(src, 'r') as zin, ZipFile(out, 'w', ZIP_DEFLATED) as zout:
    for item in zin.infolist():
        data = zin.read(item.filename)
        if item.filename == 'word/media/image24.png':
            data = replacement
        elif item.filename == 'word/document.xml':
            data = data.replace(old_text.encode('utf-8'), new_text.encode('utf-8'))
            data = data.replace(old_cap.encode('utf-8'), new_cap.encode('utf-8'))
            # Keep the replacement figure and its caption together on the page.
            start = data.find(b'<a:blip r:embed="rId32"')
            if start >= 0:
                end = data.find(b'</wp:inline>', start)
                block = data[start:end]
                block = re.sub(rb'cy="[0-9]+"', b'cy="2100000"', block)
                data = data[:start] + block + data[end:]
        zout.writestr(item, data)
print(out)
