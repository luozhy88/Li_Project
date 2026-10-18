"""Compute publication-level counts and transparent bibliometric networks."""
from pathlib import Path
import os,json,re,csv,unicodedata,collections,itertools,sys
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'output/06_完整分析';FIG=OUT/'figures';FIG.mkdir(exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/acupuncture_mpl')
import numpy as np, scipy, scipy.stats as st, networkx as nx
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':220,'svg.fonttype':'none'})
records=json.loads((OUT/'screened_records.json').read_text());rs=[r for r in records if r['decision']=='include'];N=len(rs)
def csvsave(name,rows,fields=None):
    rows=list(rows)
    with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fields or list(rows[0]));w.writeheader();w.writerows(rows)
def norm(s):
    s=unicodedata.normalize('NFKC',s).lower().strip().strip('*')
    s=re.sub('[‐‑–—−]','-',s);return re.sub(r'\s+',' ',s).strip(' .;')
def authornorm(s):
    parts=norm(s).split(',',1)
    return parts[0]+', '+parts[1].strip().replace('-','') if len(parts)==2 else parts[0]
aliases={
 'electro-acupuncture':'electroacupuncture','electro acupuncture':'electroacupuncture',
 'acupuncture therapy':'acupuncture','acupuncture treatment':'acupuncture',
 'deep dry needling':'dry needling','dry-needling':'dry needling',
 'muscle spasticity':'spasticity','muscular spasticity':'spasticity',
 'poststroke spasticity':'post-stroke spasticity','post-stroke limb spasticity':'post-stroke spasticity',
 'post-stroke muscle spasticity':'post-stroke spasticity','poststroke':'stroke','stoke':'stroke',
 'meta analysis':'meta-analysis','meta-analyses':'meta-analysis','meta-analysis.':'meta-analysis',
 'systematic reviews':'systematic review','nma':'network meta-analysis',
 'randomized controlled trial (rct)':'randomized controlled trial','randomised controlled trial':'randomized controlled trial',
 'randomized controlled trials':'randomized controlled trial','randomised controlled trials':'randomized controlled trial',
 'upper extremity':'upper limb','lower extremity':'lower limb','h reflex':'h-reflex',
 'brain plasticity':'neuroplasticity','head acupuncture':'scalp acupuncture',
 'kcc2 cotransporter':'kcc2','potassium chloride cotransporter 2':'kcc2',
 'botox':'botulinum toxin','botulinum toxin type a':'botulinum toxin',
}
keyword_map={}
for r in rs:
    r['issue_year']=int(r['year']);years=[r['issue_year']]+[int(d['value'][:4]) for d in r.get('dates',[]) if re.match(r'\d{4}',d['value'])]
    r['analysis_year']=min(years)
    terms=set()
    for raw in r['keywords']:
        parts=re.split(r'\s*[|;]\s*',raw)
        canonical=[aliases.get(norm(p),norm(p)) for p in parts if norm(p)]
        keyword_map[raw]=' | '.join(canonical);terms.update(canonical)
    r['terms']=sorted(terms)
    t=norm(r['title']);a=norm(r['abstract']);types=r['publication_types']
    if r['pmid']=='40396979':design='Consensus study'
    elif r['pmid'] in ['23869681']:design='Other human research'
    elif any(x in types for x in ['Review','Systematic Review','Meta-Analysis','Scoping Review','Network Meta-Analysis']) or re.search(r'review|meta.analysis',t) or r['pmid'] in ['33293991','10437214']:design='Review/evidence synthesis'
    elif re.search(r'\brats?\b|\brabbits?\b|\bmice\b',t):design='Animal experiment'
    elif r['pmid'] in ['41750240','15118778']:design='Secondary analysis'
    elif 'Randomized Controlled Trial' in types or re.search(r'randomi[sz]ed|randomly (?:assigned|allocated|divided)',t+' '+a):design='Randomized clinical study'
    else:design='Other human research'
    if r['pmid']=='38013386':design='Other human research' # Title calls it observational, despite contradictory abstract/indexing.
    r['design']=design
    if r['pmid']=='42264823':disease='Mixed neurological disorders'
    elif re.search('stroke|apople|cerebrovascular',t):disease='Stroke'
    elif 'cerebral palsy' in t:disease='Cerebral palsy'
    elif 'multiple sclerosis' in t:disease='Multiple sclerosis'
    elif re.search('hemifacial|facial spasm|facial paralysis|facial synkinesis',t):disease='Facial spasm/synkinesis'
    elif re.search('brain injur|disorder.*consciousness',t):disease='Brain injury/disorders of consciousness'
    elif re.search('spinal cord|htlv|myelopathy',t):disease='Spinal cord disorders'
    elif re.search('dystonia|torticollis|dysphonia',t):disease='Dystonia/spasmodic disorders'
    elif re.search('bronchi|airway|asthma|vasomotor|colonic',t):disease='Visceral/vascular spasm'
    elif re.search('pelvic floor|myofascial',t):disease='Other skeletal muscle spasm'
    else:disease='Mixed/unspecified'
    if r['pmid'] in ['40021946','10437214','1493777','37147753']:disease='Stroke'
    if r['pmid']=='40308713':disease='Other skeletal muscle spasm'
    r['disease']=disease
    r['dry_needling_focus']=bool(re.search(r'dry needl',t)) or r['pmid']=='35024159'
    r['electroacupuncture_focus']=bool(re.search(r'electro.?acup|electric acupuncture',t))
    r['author_keys']=sorted({authornorm(x) for x in r['authors']})
csvsave('04_关键词规范化词典.csv',[{'raw':k,'normalized':v} for k,v in sorted(keyword_map.items())])

# Institution-family normalization: explicit university names and their affiliates
# are grouped; hospitals without an explicit university remain separate.
uni_names=['Beijing University of Chinese Medicine','Shanghai University of Traditional Chinese Medicine','Tianjin University of Traditional Chinese Medicine','Heilongjiang University of Chinese Medicine','Fujian University of Traditional Chinese Medicine','Guangzhou University of Chinese Medicine','Changchun University of Chinese Medicine','Liaoning University of Traditional Chinese Medicine','Chengdu University of Traditional Chinese Medicine','Hunan University of Traditional Chinese Medicine','Zhejiang Chinese Medical University','Nanjing University of Chinese Medicine','Anhui University of Chinese Medicine','Henan University of Chinese Medicine','Guangxi University of Chinese Medicine','Tehran University of Medical Sciences','Iran University of Medical Sciences','Shiraz University of Medical Sciences','Ahvaz Jundishapur University of Medical Sciences','Isfahan University of Medical Sciences','Hamadan University of Medical Sciences','Arak University of Medical Sciences','Universidad Rey Juan Carlos','Universidad Complutense de Madrid','Universidad Europea de Madrid','University of Valencia','University of Cádiz','University of Sevilla','University of Zaragoza','Universidad de La Laguna','Universidad Europea de Canarias','Universitat Internacional de Catalunya','Catholic University of Ávila','Universidad Pontificia de Salamanca','Francisco de Vitoria University','San Jorge University','University of Kansas','University of Texas Health Science Center at Houston','University of Texas Medical Branch','University of Oklahoma','University of Arizona','University of Cincinnati','University of Missouri','University of South Florida','University of Maryland','Northwestern University','Kansas City University','Harvard','Cairo University','Mansoura University','Beni-Suef University','Jazan University','Majmaah University','Jouf University','Princess Nourah bint Abdulrahman University','Batterjee Medical College','Beirut Arab University','University of Chile','RMIT University','University of Melbourne','University of Antwerp','KU Leuven','University College London','University of Manchester','University of Hull','Dongguk University','Kyung Hee University','Kyung-Hee University','Kaohsiung Medical University','Chang Gung University','Chulalongkorn University','Mie University','Huazhong University of Science and Technology','Wuhan University of Science and Technology','Sichuan University','Peking University','Fuzhou University','Ningbo University','Zhengzhou University','Hebei Medical University','Chengdu Medical College','Anhui Medical University','Lanzhou University','Shenzhen University','Guangdong Pharmaceutical University','Gannan Medical University','Binzhou Medical University','Zunyi Medical University','Tongji University','City University of New York','Nanjing Medical University','Chongqing University of Traditional Chinese Medicine','Chongqing Medical University','Zhejiang University','Guizhou Medical University','Jilin University','Nanchang University','Capital Medical University','Huanggang Normal University','University Medical Center Hamburg-Eppendorf']
extra_alias={'guangzhou university of traditional chinese medicine':'Guangzhou University of Chinese Medicine','liaoning university of traditional, chinese medicine':'Liaoning University of Traditional Chinese Medicine','avila catholic university':'Catholic University of Ávila','university rey juan carlos':'Universidad Rey Juan Carlos','kyung-hee university':'Kyung Hee University','anhui college of traditional chinese medicine':'Anhui University of Chinese Medicine','hunan college of traditional chinese medicine':'Hunan University of Traditional Chinese Medicine','ucl queen square':'University College London','universidad complutense':'Universidad Complutense de Madrid','universidad de sevilla':'University of Sevilla','rey juan carlos university':'Universidad Rey Juan Carlos'}
country_patterns={'China (mainland)':r'\bchina\b|\bchn\b','Taiwan':r'\btaiwan\b','Hong Kong':r'hong kong','United States':r'\busa\b|united states|\bu\.s\.a\b','United Kingdom':r'\buk\b|united kingdom|\bengland\b','South Korea':r'\bkorea\b','Iran':r'\biran\b','Spain':r'\bspain\b','Brazil':r'\bbrazil\b','Australia':r'\baustralia\b','Japan':r'\bjapan\b','Singapore':r'\bsingapore\b','Egypt':r'\begypt\b','Saudi Arabia':r'saudi arabia','Lebanon':r'\blebanon\b','Belgium':r'\bbelgium\b','Germany':r'\bgermany\b|\bdeu\b','India':r'\bindia\b|\bind\b','Italy':r'\bitaly\b','Turkey':r'\bturkey\b|\btürkiye\b','Chile':r'\bchile\b','Thailand':r'\bthailand\b','Poland':r'\bpoland\b','Bulgaria':r'bulgarian','France':r'\bfrance\b','Portugal':r'\bportugal\b','Switzerland':r'\bswitzerland\b','Canada':r'\bcanada\b'}
affmap={}
country_patterns.update({'Norway':r'\bnorway\b','Denmark':r'\bdenmark\b','Sweden':r'\bsweden\b','Netherlands':r'\bnetherlands\b','New Zealand':r'new zealand','Malaysia':r'\bmalaysia\b','Austria':r'\baustria\b'})
extra_alias.update({'shanghai jiaotong university':'Shanghai Jiao Tong University','shanghai jiao tong university':'Shanghai Jiao Tong University','hunan university of chinese medicine':'Hunan University of Traditional Chinese Medicine','guangxi university of cm':'Guangxi University of Chinese Medicine','anhui university of cm':'Anhui University of Chinese Medicine','chinese university of hongkong':'Chinese University of Hong Kong','service hospital universitario de canarias':'Hospital Universitario de Canarias','universidad san jorge':'San Jorge University','chinese academy of sciences (cas)':'Chinese Academy of Sciences','university hospital of leuven':'KU Leuven','new york university langone':'New York University'})
for r in rs:
    r['countries']=[];r['institutions']=[]
    for aff in r['affiliations']:
        s=norm(re.sub(r'[\w.+-]+@[\w.-]+','',aff));countries={c for c,p in country_patterns.items() if re.search(p,s)}
        if 'Hong Kong' in countries or 'Taiwan' in countries:countries.discard('China (mainland)')
        orgs={u for u in uni_names if norm(u) in s}
        orgs.update(u for a,u in extra_alias.items() if a in s)
        orgs.discard('Kyung-Hee University')
        if not orgs:
            chunks=[v.strip(' .0123456789') for v in re.split('[,;]',aff) if re.search(r'hospital|academy|institute|physiocare|clinic|univers|college',v,re.I) and not re.search('e-mail|@',v)]
            chunks=[c for c in chunks if norm(c) not in {'clinician','college of chinese medicine','graduate institute of acupuncture science','university ave'}]
            if chunks:orgs.add(re.sub(r'^(?:the )','',chunks[-1],flags=re.I))
        # Location-only completions are explicit and listed in the mapping export.
        completion='explicit country/region'
        if not countries:
            if any(re.search('chinese medicine',u,re.I) for u in orgs) or re.search('beijing|changsha|ningbo|anhui',s):countries={'China (mainland)'};completion='unambiguous institution/location'
            elif re.search('texas medical branch|northwestern university|bethesda|missouri-columbia',s):countries={'United States'};completion='unambiguous institution/location'
        affmap[aff]={'affiliation':aff,'countries':' | '.join(sorted(countries)),'institutions':' | '.join(sorted(orgs)),'country_basis':completion if countries else 'unresolved'}
        r['countries'].extend(countries);r['institutions'].extend(orgs)
    r['countries']=sorted(set(r['countries']));r['institutions']=sorted(set(r['institutions']))
csvsave('05_机构与国家映射.csv',affmap.values())
# Preserve name variants and surface potentially ambiguous recurring names.
name_rows=[]
for r in rs:
    for raw in r['authors']:name_rows.append({'pmid':r['pmid'],'raw_name':raw,'normalized_name':authornorm(raw),'institution_families':'; '.join(r['institutions']),'rule':'Full given name; case/given-name hyphen normalization only; initials not expanded; residual homonym risk'})
csvsave('05b_作者名称规范化.csv',name_rows)
csvsave('06_纳入文献分类与字段覆盖.csv',[{'pmid':r['pmid'],'title':r['title'],'year':r['analysis_year'],'issue_year':r['issue_year'],'design':r['design'],'disease':r['disease'],'keywords':'; '.join(r['terms']),'countries':'; '.join(r['countries']),'institutions':'; '.join(r['institutions']),'dry_needling_focus':r['dry_needling_focus']} for r in rs])
def counts(field):return collections.Counter(v for r in rs for v in sorted(set(r[field])))
def makegraph(field,subset=rs,threshold=1):
    freq=collections.Counter(v for r in subset for v in sorted(set(r[field])));g=nx.Graph()
    g.add_nodes_from((v,{'frequency':n}) for v,n in sorted(freq.items()) if n>=threshold)
    for r in subset:
        for a,b in itertools.combinations(sorted(set(r[field]) & set(g)),2):
            if g.has_edge(a,b):g[a][b]['weight']+=1
            else:g.add_edge(a,b,weight=1)
    for a,b,d in g.edges(data=True):d['distance']=1/d['weight']
    return g,freq
def graphstats(g):
    if not g:return {}
    cs=nx.community.louvain_communities(g,weight='weight',seed=42) if g.number_of_edges() else [{x} for x in g]
    cs=sorted(cs,key=lambda c:(-sum(g.nodes[v]['frequency'] for v in c),sorted(c)))
    cmap={v:i+1 for i,c in enumerate(cs) for v in c}
    nx.set_node_attributes(g,cmap,'cluster')
    return {'nodes':len(g),'edges':g.number_of_edges(),'components':nx.number_connected_components(g),'largest_component':max(map(len,nx.connected_components(g))),'density':nx.density(g),'modularity':nx.community.modularity(g,cs,weight='weight') if g.number_of_edges() else None,'communities':[sorted(c,key=lambda v:(-g.nodes[v]['frequency'],v)) for c in cs]}
def exportgraph(field,threshold=1):
    g,freq=makegraph(field,threshold=threshold);stats=graphstats(g)
    bc=nx.betweenness_centrality(g,weight='distance',normalized=True)
    csvsave(f'network_{field}_nodes.csv',[{'label':n,'papers':g.nodes[n]['frequency'],'degree':g.degree(n),'total_link_strength':g.degree(n,weight='weight'),'betweenness':bc[n],'cluster':g.nodes[n]['cluster']} for n in sorted(g,key=lambda n:(-g.nodes[n]['frequency'],n))])
    csvsave(f'network_{field}_edges.csv',[{'source':a,'target':b,'weight':d['weight']} for a,b,d in g.edges(data=True)],['source','target','weight'])
    nx.write_graphml(g,OUT/f'network_{field}.graphml')
    return g,stats
graphs={};stats={}
for field,threshold in [('countries',1),('institutions',2),('author_keys',3),('terms',3)]:graphs[field],stats[field]=exportgraph(field,threshold)
freq=counts('terms');journals=collections.Counter(r['journal'] for r in rs);designs=collections.Counter(r['design'] for r in rs);diseases=collections.Counter(r['disease'] for r in rs)
for name,c in [('journals',journals),('designs',designs),('diseases',diseases),('countries',counts('countries')),('institutions',counts('institutions')),('authors',counts('author_keys')),('keywords',freq)]:
    csvsave(f'ranking_{name}.csv',[{'label':k,'papers':v,'percent_of_all':100*v/N} for k,v in sorted(c.items(),key=lambda z:(-z[1],z[0]))])
years=list(range(min(r['analysis_year'] for r in rs),2027));annual=collections.Counter(r['analysis_year'] for r in rs);kwannual=collections.Counter(r['analysis_year'] for r in rs if r['terms'])
csvsave('07_年度发文与关键词覆盖.csv',[{'year':y,'all':annual[y],'keyword_documents':kwannual[y],'keyword_coverage':kwannual[y]/annual[y] if annual[y] else ''} for y in years])
xy=np.array([y for y in years if y<=2025]);yy=np.array([annual[y] for y in xy]);fit=st.linregress(xy,yy)
csvsave('08_关键词逐年频次.csv',[{'keyword':k,**{str(y):sum(k in r['terms'] and r['analysis_year']==y for r in rs) for y in years}} for k in sorted(freq,key=lambda k:(-freq[k],k))])
sens=[]
for subset_name,subset in [('all',rs),('original',[r for r in rs if r['retrieval_group']=='original']),('primary_research',[r for r in rs if r['design']!='Review/evidence synthesis']),('stroke',[r for r in rs if r['disease']=='Stroke'])]:
    for th in [2,3,5]:
        g,f=makegraph('terms',subset,th);s=graphstats(g);sens.append({'subset':subset_name,'papers':len(subset),'keyword_papers':sum(bool(r['terms']) for r in subset),'threshold':th,**{k:v for k,v in s.items() if k!='communities'},'top10':'; '.join(k for k,v in f.most_common(10))})
csvsave('09_敏感性分析.csv',sens)
# Predefined calendar periods; compare proportions, not raw counts.
periods=[('1976-2009',1976,2009),('2010-2015',2010,2015),('2016-2020',2016,2020),('2021-2025',2021,2025),('2026 partial',2026,2026)]
period_rows=[]
for label,lo,hi in periods:
    sub=[r for r in rs if lo<=r['analysis_year']<=hi]
    period_rows.append({'period':label,'papers':len(sub),'keyword_papers':sum(bool(r['terms']) for r in sub),'animal':sum(r['design']=='Animal experiment' for r in sub),'reviews':sum(r['design']=='Review/evidence synthesis' for r in sub),'dry_needling_focus':sum(r['dry_needling_focus'] for r in sub)})
csvsave('10_时期研究构成.csv',period_rows)
before=[r for r in rs if 2016<=r['analysis_year']<=2020];after=[r for r in rs if 2021<=r['analysis_year']<=2025]
comparisons=[]
for category,test in [('Animal experiment',lambda r:r['design']=='Animal experiment'),('Review/evidence synthesis',lambda r:r['design']=='Review/evidence synthesis'),('Dry needling focus',lambda r:r['dry_needling_focus'])]:
    n0=sum(test(r) for r in before);n1=sum(test(r) for r in after);odds,p=st.fisher_exact([[n1,len(after)-n1],[n0,len(before)-n0]])
    comparisons.append({'category':category,'earlier_n':n0,'earlier_total':len(before),'later_n':n1,'later_total':len(after),'odds_ratio':odds,'p_exploratory':p})
ordered=sorted(comparisons,key=lambda r:r['p_exploratory']);q=1.0
for i in range(len(ordered)-1,-1,-1):
    q=min(q,ordered[i]['p_exploratory']*len(ordered)/(i+1));ordered[i]['p_bh_adjusted']=q
csvsave('11_时期比例比较.csv',comparisons)
# Two-state binomial frequency-increase model; this is NOT a CiteSpace run.
bursts=[];time=list(range(min(years),2026));T=len(time);den=np.array([kwannual[y] for y in time])
for term,count in freq.items():
    num=np.array([sum(term in r['terms'] and r['analysis_year']==y for r in rs) for y in time]);p0=num.sum()/max(1,den.sum())
    if num.sum()<5 or p0<=0:continue
    probs=[p0,min(2*p0,.999)];cost=np.array([-(num*np.log(p)+(den-num)*np.log(1-p)) for p in probs]).T
    dp=np.full((T,2),np.inf);parent=np.zeros((T,2),int);pen=np.log(T)
    dp[0]=cost[0]+[0,pen]
    for i in range(1,T):
        for state in range(2):
            opts=[dp[i-1,prev]+(pen if state>prev else 0) for prev in range(2)]
            parent[i,state]=int(np.argmin(opts));dp[i,state]=min(opts)+cost[i,state]
    states=[int(np.argmin(dp[-1]))]
    for i in range(T-1,0,-1):states.append(int(parent[i,states[-1]]))
    states=states[::-1]
    for state,group in itertools.groupby(enumerate(states),lambda z:z[1]):
        ix=[i for i,s in group]
        if state and len(ix)>=2:bursts.append({'keyword':term,'start':time[ix[0]],'end':time[ix[-1]],'length':len(ix),'log_likelihood_gain':float(sum(cost[i,0]-cost[i,1] for i in ix)),'end_censored':ix[-1]==T-1})
csvsave('12_两状态词频增强区段.csv',bursts,['keyword','start','end','length','log_likelihood_gain','end_censored'])
summary={'included':N,'excluded':len(records)-N,'original_included':sum(r['retrieval_group']=='original' for r in rs),'additional_included':sum(r['retrieval_group']=='additional' for r in rs),'years':[min(years),max(years)],'keyword_papers':sum(bool(r['terms']) for r in rs),'affiliation_papers':sum(bool(r['affiliations']) for r in rs),'country_papers':sum(bool(r['countries']) for r in rs),'institution_papers':sum(bool(r['institutions']) for r in rs),'unique_authors':len(counts('author_keys')),'unique_keywords':len(freq),'unique_journals':len(journals),'unique_countries':len(counts('countries')),'unique_institutions':len(counts('institutions')),'international_papers':sum(len(r['countries'])>1 for r in rs),'network_stats':stats,'designs':dict(designs),'diseases':dict(diseases),'periods':period_rows,'annual':dict(sorted(annual.items())),'top_keywords':freq.most_common(20),'top_journals':journals.most_common(10),'top_countries':counts('countries').most_common(10),'top_institutions':counts('institutions').most_common(10),'top_authors':counts('author_keys').most_common(10),'regression':{'start':min(xy).item(),'end':2025,'slope':fit.slope,'intercept':fit.intercept,'r2':fit.rvalue**2,'note':'Descriptive OLS; all calendar years including zeros; partial 2026 excluded'},'references_papers':sum(bool(r.get('references')) for r in rs),'bursts':bursts,'software':{'python':sys.version.split()[0],'networkx':nx.__version__,'numpy':np.__version__,'scipy':scipy.__version__,'matplotlib':matplotlib.__version__},'exclusion_counts':dict(collections.Counter(r['exclusion_reason'] for r in records if r['decision']=='exclude'))}
(OUT/'analysis_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2));(OUT/'analysis_records.json').write_text(json.dumps(rs,ensure_ascii=False,indent=2))
if '--no-figures' in sys.argv:
    print('Recomputed all statistics without rendering figures.');sys.exit(0)
def save(fig,name):
    fig.savefig(FIG/(name+'.png'),bbox_inches='tight');fig.savefig(FIG/(name+'.svg'),bbox_inches='tight');plt.close(fig)
colors=['#255D81','#2A9D8F','#E9B44C','#C95D63','#7D6B91','#7E9D75','#D98E62']
# Flowchart: mutually exclusive reasons, not a fabricated two-reviewer PRISMA process.
fig,ax=plt.subplots(figsize=(11,7));ax.axis('off')
boxes=[(.05,.74,.38,.18,'Archived expanded PubMed search\nn = 604 (395 original + 209 additional)'),(.05,.43,.38,.18,'Record-level eligibility assessment\nn = 604; unique PMID = 604'),(.05,.08,.38,.18,f'Included publications\nn = {N}\nKeyword subset: {summary["keyword_papers"]}'),(.52,.20,.44,.65,'Excluded: '+str(604-N)+'\n\n'+'\n'.join(f'{k}: {v}' for k,v in summary['exclusion_counts'].items()))]
for x,y,w,h,t in boxes:
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.015',facecolor='#F1F6F8',edgecolor=colors[0]));ax.text(x+w/2,y+h/2,t,ha='center',va='center',fontsize=10)
for start,end in [((.24,.74),(.24,.61)),((.24,.43),(.24,.26)),((.43,.52),(.52,.52))]:ax.annotate('',xy=end,xytext=start,arrowprops={'arrowstyle':'->','color':colors[0],'lw':1.5})
save(fig,'Fig01_screening')
fig,ax=plt.subplots(figsize=(11,4.5));ax.bar(years,[annual[y] for y in years],color=[colors[0] if y<2026 else colors[2] for y in years]);ax.plot(xy,fit.intercept+fit.slope*xy,color=colors[3],lw=1.5,label=f'Descriptive linear fit, R² = {fit.rvalue**2:.3f}');ax.set(xlabel='Earliest recorded publication year',ylabel='Publications');ax.legend(frameon=False);ax.set_xlim(min(years)-1,2027);ax.text(2026,annual[2026]+.5,'Partial',ha='center',fontsize=9);save(fig,'Fig02_annual')
def barplot(counter,name,xlabel='Publications (full counting)',n=10):
    vals=sorted(counter.items(),key=lambda kv:(-kv[1],kv[0]))[:n][::-1];fig,ax=plt.subplots(figsize=(10,5));ax.barh([k for k,v in vals],[v for k,v in vals],color=colors[0]);ax.set_xlabel(xlabel)
    for i,(k,v) in enumerate(vals):ax.text(v+.15,i,str(v),va='center')
    ax.set_xlim(0,max(v for k,v in vals)*1.14);save(fig,name)
barplot(journals,'Fig03_journals')
def networkplot(g,name,title,maxlabels=25):
    fig,ax=plt.subplots(figsize=(12,8));pos=nx.spring_layout(g,seed=42,k=2.5/max(1,np.sqrt(len(g))),iterations=250,weight=None)
    nx.draw_networkx_edges(g,pos,ax=ax,alpha=.25,width=[.35+np.log1p(d['weight'])*.6 for a,b,d in g.edges(data=True)],edge_color='#607D8B')
    nx.draw_networkx_nodes(g,pos,ax=ax,node_size=[40+g.nodes[n]['frequency']*30 for n in g],node_color=[colors[(g.nodes[n]['cluster']-1)%len(colors)] for n in g],alpha=.9,linewidths=.7,edgecolors='white')
    labels={n:n for n in sorted(g,key=lambda n:(-g.nodes[n]['frequency'],n))[:maxlabels]}
    # Place labels in two ordered outside columns to avoid unreadable overlaps.
    for side in [-1,1]:
        nodes=sorted([n for n in labels if (pos[n][0]<0)==(side<0)],key=lambda n:pos[n][1])
        for n,y in zip(nodes,np.linspace(-1.08,1.08,len(nodes))):
            ax.annotate(n,xy=pos[n],xytext=(side*1.45,y),fontsize=8,ha='right' if side<0 else 'left',va='center',arrowprops={'arrowstyle':'-','color':'#9AA9B0','lw':.5},bbox={'facecolor':'white','edgecolor':'none','alpha':.85,'pad':1})
    ax.set_xlim(-2.5,2.5);ax.set_ylim(-1.25,1.25);ax.set_title(title,pad=15);ax.axis('off');save(fig,name)
networkplot(graphs['countries'],'Fig04_countries','Country/region coauthorship | all nodes',30)
barplot(counts('institutions'),'Fig05_institutions',n=10)
networkplot(graphs['author_keys'],'Fig06_authors','Coauthor name network | at least 3 publications',24)
networkplot(graphs['terms'],'Fig07_keywords','Keyword co-occurrence | at least 3 publications',45)
top=[k for k,v in freq.most_common(20)];heatyears=list(range(2010,2027));data=np.array([[sum(k in r['terms'] and r['analysis_year']==y for r in rs)/kwannual[y]*100 if kwannual[y] else np.nan for y in heatyears] for k in top]);fig,ax=plt.subplots(figsize=(11,8));im=ax.imshow(np.ma.masked_invalid(data),aspect='auto',cmap='YlGnBu',vmin=0);ax.set_yticks(range(len(top)),top);ax.set_xticks(range(len(heatyears)),[str(y) if y<2026 else '2026*' for y in heatyears],rotation=60);fig.colorbar(im,ax=ax,label='% of keyword-bearing publications in the year',shrink=.7);ax.set_title('Normalized annual keyword prevalence | * partial year');save(fig,'Fig08_keyword_heatmap')
fig,axes=plt.subplots(1,2,figsize=(13,5),gridspec_kw={'width_ratios':[1,1.5]})
for ax,c,title in zip(axes,[designs,diseases],['Study designs','Disease domains']):
    items=c.most_common()[::-1];ax.barh([k for k,v in items],[v for k,v in items],color=colors[1]);ax.set_xlabel('Publications');ax.set_title(title)
    for i,(k,v) in enumerate(items):ax.text(v+.3,i,str(v),va='center',fontsize=9)
    ax.set_xlim(0,max(v for k,v in items)*1.15)
fig.tight_layout(w_pad=3);save(fig,'Fig09_design_disease')
fig,axes=plt.subplots(1,2,figsize=(11,4.5));x=np.arange(len(period_rows));width=.24
for i,(key,label) in enumerate([('animal','Animal experiments'),('reviews','Reviews/evidence synthesis'),('dry_needling_focus','Dry needling in title')]):axes[0].bar(x+(i-1)*width,[100*p[key]/p['papers'] for p in period_rows],width,color=colors[i],label=label)
axes[0].set_xticks(x,[p['period'] for p in period_rows],rotation=30,ha='right');axes[0].set_ylabel('% of period publications');axes[0].legend(fontsize=8,frameon=False);axes[0].set_title('Study composition (categories overlap)')
orig=[r for r in rs if r['retrieval_group']=='original'];add=[r for r in rs if r['retrieval_group']=='additional'];axes[1].bar(['Original eligible','Added eligible'],[len(orig),len(add)],color=[colors[0],colors[1]]);axes[1].set_ylabel('Publications');axes[1].set_title('Search-strategy contribution')
for i,g in enumerate([orig,add]):axes[1].text(i,len(g)+1,str(len(g)),ha='center')
fig.tight_layout();save(fig,'Fig10_period_sensitivity')
print(json.dumps(summary,ensure_ascii=False,indent=2))
