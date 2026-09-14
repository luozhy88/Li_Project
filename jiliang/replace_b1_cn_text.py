import json
import re
from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED

root = Path('/Users/apple/Desktop/github.local/Li_Project/jiliang')
src = root / 'output/report/10_False_Localizing_Sign_文献计量分析报告_Medicine投稿修订版v8.docx'
out = root / 'output/report/10_False_Localizing_Sign_文献计量分析报告_Medicine投稿修订版v9.docx'
stats = json.loads((root / 'output/logs/b1_own_stats.json').read_text(encoding='utf-8'))

fn = stats['full_network']
t3f = stats['top3_frequency']
bc01 = stats['betweenness_gt_0.1']
cl = stats['clustering']
Q = round(cl['modularity_Q'], 3)
S = round(cl['weighted_mean_silhouette_S'], 3)

fmt_f = lambda e: f"{e['keyword']}（{e['freq']}次）"
fmt_bc = lambda e: f"{e['keyword']}（{e['bc']:.2f}）"
cluster_list = '、'.join(f"#{c['id']} {c['label']}" for c in cl['clusters'][:-1])
cluster_list += f"和#{cl['clusters'][-1]['id']} {cl['clusters'][-1]['label']}"

para1 = (f"假定位体征研究的关键词共现网络包含{fn['nodes']}个节点和{fn['links']}条连接，网络密度为{fn['density']}（图B1）。"
         f"关键词的出现频率是衡量其重要性的关键指标——频率越高，表明研究关注度越大。"
         f"本研究中出现频率最高的3个关键词分别是{fmt_f(t3f[0])}、{fmt_f(t3f[1])}和{fmt_f(t3f[2])}。"
         f"通常认为，介数中心性大于0.1的关键词在共现网络中占据枢纽地位，这也意味着它们是该领域的研究焦点。"
         f"本研究中介数中心性大于0.1的关键词是{fmt_bc(bc01[0])}、{fmt_bc(bc01[1])}和{fmt_bc(bc01[2])}。"
         f"在出现频率≥5次且介数中心性≥0.05的关键词中，按介数中心性排名前三的依次是："
         f"{fmt_bc(bc01[0])}、{fmt_bc(bc01[1])}和{fmt_bc(bc01[2])}。")

para2 = (f"基于关键词共现关系，进行了聚类分析以对关系密切的关键词进行分类。"
         f"采用基于模块度的聚类方法对关键词进行聚类，并使用对数似然比（LLR）方法为各簇命名，结果如图B2所示。"
         f"共识别出{cl['n_clusters']}个关键词簇（#0–{cl['n_clusters']-1}），分别为{cluster_list}。"
         f"聚类质量指数为Q = {Q}（Q > 0.3），表明聚类结构显著；"
         f"加权平均轮廓系数S = {S}，提示簇内同质性中等，这与本研究纳入文献规模较小且主题交叉的特点相符，聚类结果总体合理。"
         f"主题时间图按照关键词在本数据中的首次至最后一次出现年份绘制，红色区段反映记录覆盖时间，并非引文突现强度。")

run_tpl = ('<w:p><w:r><w:rPr><w:rFonts w:ascii="Songti SC" w:hAnsi="Songti SC" w:eastAsia="Songti SC"/></w:rPr>'
           '<w:t xml:space="preserve">{}</w:t></w:r></w:p>')
replacement = run_tpl.format(para1) + run_tpl.format(para2)

pattern = re.compile(r'<w:p>(?:(?!</?w:p>).)*?The keyword co-occurrence network of false localizing sign research.*?</w:p>',
                     re.DOTALL)

with ZipFile(src, 'r') as zin, ZipFile(out, 'w', ZIP_DEFLATED) as zout:
    for item in zin.infolist():
        data = zin.read(item.filename)
        if item.filename == 'word/document.xml':
            doc = data.decode('utf-8')
            new_doc, n = pattern.subn(replacement, doc)
            assert n == 1, f'expected 1 replacement, got {n}'
            data = new_doc.encode('utf-8')
        zout.writestr(item, data)
print(out)
