from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

root = Path('/Users/apple/Desktop/github.local/Li_Project/jiliang')
src = root / 'output/report/10_False_Localizing_Sign_文献计量分析报告_Medicine投稿修订版v4.docx'
out = root / 'output/report/10_False_Localizing_Sign_文献计量分析报告_Medicine投稿修订版v6.docx'
fig = root / 'input/info/false_localizing_info_style_figures/fig10_keyword_network_and_clusters.png'

old_text = ('关键词网络显示影像学与鉴别诊断、颅内压相关表现、脑肿瘤/占位性病变、颈髓/脊髓压迫以及颅神经和运动体征等相互关联的主题。主题时间图按照关键词在本数据中的首次至最后一次出现年份绘制，红色区段反映记录覆盖时间，并非引文突现强度。')
new_text = ('Using the 159 records included in this study, the keyword co-occurrence network was reconstructed from MeSH terms and other keywords appearing in the same record. For readability, the visualization retains the 26 most frequent substantive terms and 74 co-occurrence links with link weight >= 2; node size is proportional to keyword frequency. The resulting communities highlight imaging and differential diagnosis, intracranial pressure and CSF-related conditions, brain or meningeal lesions, spinal or cervical localization, and cranial nerve or motor signs. Panel A shows the frequency-scaled co-occurrence network, whereas Panel B colors the same data by greedy-modularity communities. These are descriptive network communities, not citation-based clusters.')
old_cap = 'Fig. B1. Keyword co-occurrence network (A) and modularity-based clustered thematic structure (B). Synonymous variants were merged using a thesaurus file before analysis.'
new_cap = ('Fig. B1. Keyword co-occurrence network (A) and data-driven community visualization (B) based on the 159 records included in this study. '
            'The graph displays 26 substantive keywords and 74 retained links (co-occurrence weight >= 2); node size represents frequency and colors in panel B represent greedy-modularity communities.')

with ZipFile(src, 'r') as zin, ZipFile(out, 'w', ZIP_DEFLATED) as zout:
    fig_bytes = fig.read_bytes()
    for item in zin.infolist():
        data = zin.read(item.filename)
        if item.filename == 'word/media/image24.png':
            data = fig_bytes
        elif item.filename == 'word/document.xml':
            data = data.replace(old_text.encode(), new_text.encode())
            data = data.replace(old_cap.encode(), new_cap.encode())
        zout.writestr(item, data)
print(out)
