from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import json

from PIL import Image

root = Path('/Users/apple/Desktop/github.local/Li_Project/jiliang')
src = root / 'output/report/10_False_Localizing_Sign_文献计量分析报告_Medicine投稿修订版v4.docx'
out = root / 'output/report/10_False_Localizing_Sign_文献计量分析报告_Medicine投稿修订版v8.docx'
fig_net = root / 'output/figures/figB1_own_keyword_network.png'
fig_clu = root / 'output/figures/figB2_own_keyword_clusters.png'
stats = json.loads((root / 'output/logs/b1_own_stats.json').read_text(encoding='utf-8'))

fn = stats['full_network']
pv = stats['pruned_view']
t3f = stats['top3_frequency']
bc01 = stats['betweenness_gt_0.1']
t3b = stats['top3_betweenness_freq_ge5_bc_ge0.05']
cl = stats['clustering']
cl['weighted_mean_silhouette_S'] = round(cl['weighted_mean_silhouette_S'], 3)
cl['modularity_Q'] = round(cl['modularity_Q'], 3)

def cx_cy(img_path, cx=5715000):
    with Image.open(img_path) as im:
        w, h = im.size
    return cx, round(cx * h / w)

cx1, cy1 = cx_cy(fig_net)
cx2, cy2 = cx_cy(fig_clu)

fmt_bc = lambda e: f"{e['keyword']} ({e['bc']:.2f})"
fmt_f = lambda e: f"{e['keyword']} ({e['freq']})"
cluster_list = '; '.join(f"#{c['id']} {c['label']}" for c in cl['clusters'])

old_text = ('关键词网络显示影像学与鉴别诊断、颅内压相关表现、脑肿瘤/占位性病变、颈髓/脊髓压迫以及颅神经和运动体征等相互关联的主题。主题时间图按照关键词在本数据中的首次至最后一次出现年份绘制，红色区段反映记录覆盖时间，并非引文突现强度。')
new_text = (f"The keyword co-occurrence network of false localizing sign research was reconstructed from the {stats['records']} included records "
            f"and contains {fn['nodes']} nodes and {fn['links']:,} links, with a network density of {fn['density']} (Fig. B1). "
            f"Keyword frequency reflects the level of research attention: the three most frequent keywords were "
            f"{fmt_f(t3f[0])}, {fmt_f(t3f[1])}, and {fmt_f(t3f[2])}. "
            f"Keywords with betweenness centrality greater than 0.1 occupy hub positions in the co-occurrence network and represent "
            f"research foci; these were {fmt_bc(bc01[0])}, {fmt_bc(bc01[1])}, and {fmt_bc(bc01[2])}. "
            f"Among keywords with frequency >= 5 and betweenness centrality >= 0.05, the top three ranked by centrality were "
            f"{fmt_bc(t3b[0])}, {fmt_bc(t3b[1])}, and {fmt_bc(t3b[2])}. "
            f"Based on these co-occurrence relationships, modularity-based clustering identified {cl['n_clusters']} keyword clusters, "
            f"labeled by the log-likelihood ratio (LLR) method (Fig. B2): {cluster_list}. "
            f"The modularity (Q = {cl['modularity_Q']} > 0.3) indicates a significant cluster structure, and the weighted mean "
            f"silhouette (S = {cl['weighted_mean_silhouette_S']}) suggests moderate within-cluster homogeneity, which is expected "
            f"given the small and cross-disciplinary corpus. "
            f"The topic time map is drawn from the first to the last appearance year of each keyword in this dataset; "
            f"red segments reflect record coverage time, not citation burst strength.")

old_cap = 'Fig. B1. Keyword co-occurrence network (A) and modularity-based clustered thematic structure (B). Synonymous variants were merged using a thesaurus file before analysis.'
new_cap = (f"Fig. B1. Keyword co-occurrence network of false localizing sign research ({fn['nodes']} nodes, {fn['links']:,} links; "
           f"density = {fn['density']}). The visualization retains keywords with frequency >= {pv['min_freq']} and links with "
           f"co-occurrence weight >= {pv['min_link_weight']}; node size is proportional to keyword frequency, and purple rings "
           f"mark keywords with betweenness centrality > 0.1. Synonymous variants were merged using a thesaurus file before analysis.")
new_cap2 = (f"Fig. B2. Modularity-based cluster map of keywords in false localizing sign research. "
            f"{cl['n_clusters']} clusters (#0-{cl['n_clusters']-1}) were identified and labeled by the log-likelihood ratio (LLR) "
            f"method (modularity Q = {cl['modularity_Q']}, weighted mean silhouette S = {cl['weighted_mean_silhouette_S']}).")

new_fig2_xml = (
    '<w:p><w:pPr><w:jc w:val="center"/></w:pPr><w:r><w:drawing>'
    '<wp:inline xmlns:a="http://schemas.openxmlformats.org/drawingml/2006/main" xmlns:pic="http://schemas.openxmlformats.org/drawingml/2006/picture">'
    f'<wp:extent cx="{cx2}" cy="{cy2}"/><wp:docPr id="29" name="Picture 29"/>'
    '<wp:cNvGraphicFramePr><a:graphicFrameLocks noChangeAspect="1"/></wp:cNvGraphicFramePr>'
    '<a:graphic><a:graphicData uri="http://schemas.openxmlformats.org/drawingml/2006/picture">'
    '<pic:pic><pic:nvPicPr><pic:cNvPr id="0" name="figB2_own_keyword_clusters.png"/><pic:cNvPicPr/></pic:nvPicPr>'
    '<pic:blipFill><a:blip r:embed="rId36"/><a:stretch><a:fillRect/></a:stretch></pic:blipFill>'
    f'<pic:spPr><a:xfrm><a:off x="0" y="0"/><a:ext cx="{cx2}" cy="{cy2}"/></a:xfrm>'
    '<a:prstGeom prst="rect"/></pic:spPr></pic:pic></a:graphicData></a:graphic></wp:inline>'
    '</w:drawing></w:r></w:p>'
    '<w:p><w:pPr><w:jc w:val="center"/></w:pPr><w:r><w:rPr>'
    '<w:rFonts w:ascii="Songti SC" w:hAnsi="Songti SC" w:eastAsia="Songti SC"/><w:i/><w:sz w:val="18"/></w:rPr>'
    '<w:t>' + new_cap2 + '</w:t></w:r></w:p>')

with ZipFile(src, 'r') as zin, ZipFile(out, 'w', ZIP_DEFLATED) as zout:
    for item in zin.infolist():
        data = zin.read(item.filename)
        if item.filename == 'word/media/image24.png':
            data = fig_net.read_bytes()
        elif item.filename == 'word/document.xml':
            doc = data.decode('utf-8')
            assert old_text in doc, 'body paragraph not found'
            assert old_cap in doc, 'caption not found'
            doc = doc.replace(old_text, new_text)
            doc = doc.replace(old_cap, new_cap)
            doc = doc.replace('<wp:extent cx="5715000" cy="3554197"/><wp:docPr id="26" name="Picture 26"/>',
                              f'<wp:extent cx="{cx1}" cy="{cy1}"/><wp:docPr id="26" name="Picture 26"/>')
            doc = doc.replace('<a:ext cx="5715000" cy="3554197"/>', f'<a:ext cx="{cx1}" cy="{cy1}"/>', 1)
            doc = doc.replace('name="fig10_keyword_network_and_clusters.png"',
                              'name="figB1_own_keyword_network.png"')
            # renumber following B figures first (B3 before B2 to avoid collision)
            assert doc.count('Fig. B3.') == 1 and doc.count('Fig. B2.') == 1
            doc = doc.replace('Fig. B3.', 'Fig. B4.')
            doc = doc.replace('Fig. B2.', 'Fig. B3.')
            anchor = '<w:t>' + new_cap + '</w:t></w:r></w:p>'
            assert anchor in doc, 'B1 caption anchor not found'
            doc = doc.replace(anchor, anchor + new_fig2_xml)
            data = doc.encode('utf-8')
        elif item.filename == 'word/_rels/document.xml.rels':
            rels = data.decode('utf-8')
            rel = ('<Relationship Id="rId36" '
                   'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/image" '
                   'Target="media/image36.png"/>')
            rels = rels.replace('</Relationships>', rel + '</Relationships>')
            data = rels.encode('utf-8')
        zout.writestr(item, data)
    zout.writestr('word/media/image36.png', fig_clu.read_bytes())
print(out)
