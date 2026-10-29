"""Audit existing results and redraw descriptive figures; no enrichment or docking rerun.
Run from the project root. Outputs are confined to this dated support directory.
"""
from pathlib import Path
import hashlib
import json
import textwrap
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
B = ROOT / 'v20241015'
GROUPS = ['ZBM', 'XKC', 'XYT']
COLORS = ['#ce6548', '#208e89', '#7161a8']
sources = [ROOT/'edit.Paper/output/article_jcpsp.docx']
disease_file = B/'input/GeneCards-SearchResults_Thyroid_tumor.csv'
sources.append(disease_file)
disease = set(pd.read_csv(disease_file)['Gene Symbol'].dropna())
counts, herbs, terms = [], [], []
for g in GROUPS:
    f = B/f'output/01_TCM_CID_compound_Gene_{g}.csv'
    sources.append(f)
    df = pd.read_csv(f)
    genes = set(';'.join(df.entrez_gene_symbol.dropna()).split(';'))
    overlap = genes & disease
    common = B/f'output/02_common_gene_{g}.csv'
    sources.append(common)
    assert overlap == set(pd.read_csv(common).Com_Gene)
    counts.append(dict(group=g, disease_total=len(disease), ingredient_gene_total=len(genes),
                       overlap=len(overlap), disease_only=len(disease-genes), ingredient_only=len(genes-disease),
                       overlap_pct_of_ingredient_genes=100*len(overlap)/len(genes),
                       overlap_pct_of_disease=100*len(overlap)/len(disease),
                       unique_CID=df.CID.nunique()))
    for herb, rows in df.groupby('herb_cn_name'):
        herbs.append(dict(group=g, herb=herb, unique_CID=rows.CID.nunique(), rows=len(rows)))
    for kind in ['go','kegg']:
        f = B/f'output/03_{kind}_enrichment_results_{g}.csv'
        sources.append(f)
        e = pd.read_csv(f)
        assert (e['p.adjust'] < .05).all()
        e['group'] = g
        e['analysis'] = kind
        e['ratio_numeric'] = e.GeneRatio.map(lambda x: int(x.split('/')[0])/int(x.split('/')[1]))
        terms.append(e)
counts = pd.DataFrame(counts)
counts.to_csv(OUT/'20261023_target_counts.csv',index=False)
pd.DataFrame(herbs).to_csv(OUT/'20261023_herb_inventory.csv',index=False)
allterms = pd.concat(terms,ignore_index=True)
allterms.to_csv(OUT/'20261023_existing_enrichment_results.csv',index=False)

# Figure 1: explicit disjoint Venn-region counts, avoiding rounded union percentages.
fig, axes = plt.subplots(1,3,figsize=(10.5,3.5))
from matplotlib.patches import Circle
for ax, row, color in zip(axes, counts.to_dict('records'), COLORS):
    ax.add_patch(Circle((.39,.52),.29,facecolor='#cfdce5',alpha=.7))
    ax.add_patch(Circle((.68,.52),.23,facecolor=color,alpha=.4))
    ax.text(.25,.52,f"{row['disease_only']:,}",ha='center')
    ax.text(.52,.52,f"{row['overlap']:,}",ha='center',weight='bold')
    ax.text(.76,.52,f"{row['ingredient_only']:,}",ha='center')
    ax.text(.28,.16,'Disease only',ha='center',fontsize=9)
    ax.text(.77,.16,'Ingredient only',ha='center',fontsize=9)
    ax.set_title(row['group'] + ('*' if row['group']=='XYT' else ''))
    ax.set(xlim=(0,1),ylim=(0,1),aspect='equal'); ax.axis('off')
fig.text(.5,.02,'Disease gene set: 12,987 genes. Circle areas are not proportional to set sizes.  *Provisional dataset label.',ha='center',fontsize=9)
fig.tight_layout(rect=(0,.05,1,1))
fig.savefig(OUT/'20261023_Figure1.png',dpi=300,facecolor='white');plt.close(fig)

# Figure 2: union of three lowest adjusted-P GO terms per group and ontology.
go = allterms[allterms.analysis=='go']
fig,axes=plt.subplots(3,1,figsize=(10,11.5))
selected=[]
for ax,ont in zip(axes,['BP','MF','CC']):
    sub=go[go.ONTOLOGY==ont]
    ids=list(dict.fromkeys(sub.sort_values('p.adjust').groupby('group',sort=False).head(3).sort_values('p.adjust').ID))
    chosen=sub[sub.ID.isin(ids)].copy();selected.append(chosen)
    names=[sub[sub.ID==i].iloc[0].Description for i in ids]
    y=np.arange(len(ids))
    for j,(g,c) in enumerate(zip(GROUPS,COLORS)):
        vals=chosen[chosen.group==g].set_index('ID').ratio_numeric.reindex(ids)
        ax.barh(y+(j-1)*.24,vals,height=.22,label=g+('*' if g=='XYT' else ''),color=c)
    ax.set_yticks(y,['\n'.join(textwrap.wrap(n,48)) for n in names],fontsize=11)
    ax.invert_yaxis();ax.set_xlabel('GeneRatio');ax.set_title(ont);ax.grid(axis='x',alpha=.2)
    ax.spines[['top','right']].set_visible(False)
axes[-1].legend(loc='lower right',fontsize=8)
fig.tight_layout();fig.savefig(OUT/'20261023_Figure2.png',dpi=300,facecolor='white');plt.close(fig)
pd.concat(selected).to_csv(OUT/'20261023_Figure2_source.csv',index=False)

# Figure 3: union of five lowest adjusted-P pathways per group, plus hsa05216.
kegg=allterms[allterms.analysis=='kegg']
ids=list(dict.fromkeys(kegg.sort_values('p.adjust').groupby('group',sort=False).head(5).sort_values('p.adjust').ID))
if 'hsa05216' not in ids: ids.append('hsa05216')
chosen=kegg[kegg.ID.isin(ids)].copy()
fig,ax=plt.subplots(figsize=(10,6.5))
for j,(g,c) in enumerate(zip(GROUPS,COLORS)):
    vals=chosen[chosen.group==g].set_index('ID').ratio_numeric.reindex(ids)
    ax.barh(np.arange(len(ids))+(j-1)*.24,vals,height=.22,label=g+('*' if g=='XYT' else ''),color=c)
names=[kegg[kegg.ID==i].iloc[0].Description for i in ids]
ax.set_yticks(range(len(ids)),['\n'.join(textwrap.wrap(n,46)) for n in names],fontsize=9)
ax.invert_yaxis();ax.set_xlabel('GeneRatio');ax.legend();ax.grid(axis='x',alpha=.2)
ax.spines[['top','right']].set_visible(False)
fig.tight_layout();fig.savefig(OUT/'20261023_Figure3.png',dpi=300,facecolor='white');plt.close(fig)
chosen.to_csv(OUT/'20261023_Figure3_source.csv',index=False)

dock=[]
for folder in sorted((ROOT/'add.docking/input').glob('docking_*')):
    if not folder.is_dir():continue
    f=folder/'docking_scores.csv'; sources.append(f)
    df=pd.read_csv(f)
    code=folder.name.split('_')[1]
    het=folder/'receptor/ref_ligand.pdb';sources.append(het)
    ligands=sorted({l[17:20].strip() for l in het.read_text().splitlines() if l.startswith(('HETATM','ATOM'))})
    log=folder/'smina.log';sources.append(log)
    dock.append(dict(PDB=code,engine='smina (archived log)',poses=len(df),best_score=df.minimizedAffinity.min(),
                     autobox_reference_residue_names=';'.join(ligands),
                     status={'4P5A':'EXCLUDED: bacterial ThyX/PolB; not human POLB',
                             '6UNI':'NOT VALIDATED: heme removed; autobox includes HEM and QDV',
                             '9AYG':'NOT VALIDATED: apo structure; mixed non-protein autobox reference'}[code]))
pd.DataFrame(dock).to_csv(OUT/'20261023_docking_audit.csv',index=False)
for f in ['01_Get_compound.R','02_venny_gene.R','03_KEGG.R','03.2_KEGG_Plot.R','04_download_PDB_by_gene.R','07.venny_vina.R']:
    sources.append(B/f)
manifest=[dict(path=str(p.relative_to(ROOT)),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size) for p in sources]
(OUT/'20261023_source_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
print(counts.to_string(index=False))
print('Outputs:',OUT)
