"""
针灸穴位-基因关联数据库 Web分析平台 v2.1
基于 Python http.server + jinja2 + pandas(xlsx数据源)

启动方式:
    python src/web/app.py
    然后访问 http://localhost:8080
"""
import os
import sys
import json
import urllib.parse
from http.server import HTTPServer, BaseHTTPRequestHandler
from jinja2 import Environment, FileSystemLoader

import pandas as pd

# 路径设置
BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
TEMPLATE_DIR = os.path.join(os.path.dirname(__file__), "templates")
XLSX_PATH = os.path.join(BASE_DIR, "output", "acupoint_gene_database.xlsx")
sys.path.insert(0, os.path.join(BASE_DIR, "src"))

from web.knowledge_base import get_gene_function, get_acupoint_info

# Jinja2 环境
os.makedirs(TEMPLATE_DIR, exist_ok=True)
jenv = Environment(loader=FileSystemLoader(TEMPLATE_DIR))

# ========================================================================
# 加载 xlsx 数据到内存（应用启动时一次性加载）
# ========================================================================

print("📂 正在加载 Excel 数据源...")
_xl = pd.ExcelFile(XLSX_PATH)

# 核心数据表
df_articles = pd.read_excel(_xl, sheet_name="文献主表")
df_genes = pd.read_excel(_xl, sheet_name="基因统计")
df_acupoints = pd.read_excel(_xl, sheet_name="穴位统计")
df_diseases = pd.read_excel(_xl, sheet_name="疾病统计")
df_ag = pd.read_excel(_xl, sheet_name="文献-基因关联")      # article-gene
df_aa = pd.read_excel(_xl, sheet_name="文献-穴位关联")      # article-acupoint
df_ad = pd.read_excel(_xl, sheet_name="文献-疾病关联")      # article-disease
df_network = pd.read_excel(_xl, sheet_name="共现网络")       # 预计算共现

# 预计算常用索引，加速查询
gene_pmids_map = df_ag.groupby("基因符号")["PMID"].apply(lambda s: set(s.unique())).to_dict()
acupoint_pmids_map = df_aa.groupby("穴位英文名")["PMID"].apply(lambda s: set(s.unique())).to_dict()
disease_pmids_map = df_ad.groupby("疾病")["PMID"].apply(lambda s: set(s.unique())).to_dict()

# PMID -> 文章信息字典
pmid_to_article = {}
for _, row in df_articles.iterrows():
    pmid_to_article[row["PMID"]] = {
        "pmid": row["PMID"],
        "title": row["标题"],
        "year": int(row["年份"]) if pd.notna(row["年份"]) else None,
        "doi": row["DOI"] if pd.notna(row["DOI"]) else "",
        "journal": row["期刊"] if pd.notna(row["期刊"]) else "",
        "authors": row["作者"] if pd.notna(row["作者"]) else "",
        "abstract": row["摘要"] if pd.notna(row["摘要"]) else "",
        "has_fulltext": bool(row["有全文"]) if pd.notna(row["有全文"]) else False,
    }

# ID -> PMID 映射
id_to_pmid = dict(zip(df_articles["ID"], df_articles["PMID"]))
pmid_to_id = {v: k for k, v in id_to_pmid.items()}

print(f"   文献: {len(df_articles)} 篇")
print(f"   基因: {len(df_genes)} 个")
print(f"   穴位: {df_aa['穴位英文名'].nunique()} 个唯一名 / {len(df_aa)} 条关联")
print(f"   疾病: {len(df_diseases)} 个")
print(f"   预计算共现对: {len(df_network)} 对")
print("✅ Excel 数据加载完成\n")


# ========================================================================
# 数据查询函数（基于 pandas DataFrame）
# ========================================================================

def get_stats():
    """统计概览"""
    return {
        "articles": len(df_articles),
        "genes": len(df_genes),
        "acupoints": len(df_acupoints),
        "diseases": len(df_diseases),
        "fulltext": int(df_articles["有全文"].sum()),
        "article_gene": len(df_ag),
        "article_acupoint": len(df_aa),
        "article_disease": len(df_ad),
    }


def get_year_data():
    """年份分布"""
    counts = df_articles[df_articles["年份"].notna()].groupby("年份").size().reset_index(name="count")
    counts = counts.sort_values("年份")
    return [{"year": int(r["年份"]), "count": int(r["count"])} for _, r in counts.iterrows()]


def get_top_genes(limit=20):
    """Top 基因"""
    top = df_genes.nlargest(limit, "关联文献数")
    return [{"gene": r["基因符号"], "count": int(r["关联文献数"])} for _, r in top.iterrows()]


def get_top_acupoints(limit=20):
    """Top 穴位（去重后按关联文献数排序）"""
    # 穴位统计表已经按ID去重统计了，直接使用
    top = df_acupoints.nlargest(limit, "关联文献数")
    return [{
        "name_en": r["英文名"],
        "name_cn": r["中文名"] if pd.notna(r["中文名"]) else "",
        "code": r["代码"] if pd.notna(r["代码"]) else "",
        "count": int(r["关联文献数"])
    } for _, r in top.iterrows()]


def get_top_diseases(limit=15):
    """Top 疾病"""
    top = df_diseases.nlargest(limit, "关联文献数")
    return [{"disease": r["疾病名称"], "count": int(r["关联文献数"])} for _, r in top.iterrows()]


def get_acupoint_list():
    """所有穴位列表（用于下拉选择，按英文名去重）"""
    seen = set()
    result = []
    for _, r in df_acupoints.iterrows():
        name = r["英文名"]
        if name and name not in seen:
            seen.add(name)
            result.append({
                "id": int(r["ID"]),
                "name_en": name,
                "name_cn": r["中文名"] if pd.notna(r["中文名"]) else "",
                "code": r["代码"] if pd.notna(r["代码"]) else "",
                "meridian": r["经络"] if pd.notna(r["经络"]) else "",
            })
    return result


def get_gene_list():
    """所有基因列表（用于搜索建议）"""
    return [{"id": int(r["ID"]), "gene": r["基因符号"]} for _, r in df_genes.iterrows()]


def get_disease_list():
    """所有疾病列表（用于筛选）"""
    return [{"id": int(r["ID"]), "name": r["疾病名称"]} for _, r in df_diseases.iterrows()]


def query_acupoint_genes(acupoint_name: str):
    """
    查询指定穴位关联的所有基因
    返回: [{gene, cooccur, jaccard, gene_total, apt_total, article_count, articles:[{pmid,title,year}], gene_function}]
    """
    if not acupoint_name:
        return []

    apt_pmids = acupoint_pmids_map.get(acupoint_name, set())
    if not apt_pmids:
        return []

    apt_total = len(apt_pmids)

    # 找到与这些PMID共同出现的基因
    cooccur_df = df_ag[df_ag["PMID"].isin(apt_pmids)]
    gene_cooccur = cooccur_df.groupby("基因符号")["PMID"].nunique().reset_index(name="cooccur")
    gene_cooccur = gene_cooccur.sort_values("cooccur", ascending=False)

    results = []
    for _, row in gene_cooccur.iterrows():
        gene = row["基因符号"]
        cooccur = int(row["cooccur"])
        gene_pmids = gene_pmids_map.get(gene, set())
        gene_total = len(gene_pmids)
        union = apt_total + gene_total - cooccur
        jaccard = round(cooccur / union, 4) if union > 0 else 0

        # 支持文献（共同PMID对应的文章）
        common_pmids = list(apt_pmids & gene_pmids)
        common_pmids.sort(reverse=True)  # 按PMID降序（近似按年份降序）
        article_list = []
        for pmid in common_pmids[:20]:
            art = pmid_to_article.get(pmid)
            if art:
                article_list.append({
                    "pmid": art["pmid"],
                    "title": art["title"],
                    "year": art["year"],
                    "doi": art["doi"],
                    "journal": art["journal"],
                })

        results.append({
            "gene": gene,
            "cooccur": cooccur,
            "jaccard": jaccard,
            "apt_total": apt_total,
            "gene_total": gene_total,
            "article_count": len(article_list),
            "articles": article_list,
            "gene_function": get_gene_function(gene),
        })

    return results


def query_gene_acupoints(gene_symbol: str):
    """
    查询指定基因关联的所有穴位
    返回: [{acupoint, name_cn, code, cooccur, jaccard, gene_total, apt_total, article_count, articles, meridian, apt_function}]
    """
    if not gene_symbol:
        return []

    gene_pmids = gene_pmids_map.get(gene_symbol, set())
    if not gene_pmids:
        return []

    gene_total = len(gene_pmids)

    # 找到与这些PMID共同出现的穴位
    cooccur_df = df_aa[df_aa["PMID"].isin(gene_pmids)]
    apt_cooccur = cooccur_df.groupby("穴位英文名")["PMID"].nunique().reset_index(name="cooccur")
    apt_cooccur = apt_cooccur.sort_values("cooccur", ascending=False)

    results = []
    for _, row in apt_cooccur.iterrows():
        acupoint = row["穴位英文名"]
        cooccur = int(row["cooccur"])
        apt_pmids = acupoint_pmids_map.get(acupoint, set())
        apt_total = len(apt_pmids)
        union = gene_total + apt_total - cooccur
        jaccard = round(cooccur / union, 4) if union > 0 else 0

        # 支持文献
        common_pmids = list(gene_pmids & apt_pmids)
        common_pmids.sort(reverse=True)
        article_list = []
        for pmid in common_pmids[:20]:
            art = pmid_to_article.get(pmid)
            if art:
                article_list.append({
                    "pmid": art["pmid"],
                    "title": art["title"],
                    "year": art["year"],
                    "doi": art["doi"],
                    "journal": art["journal"],
                })

        # 穴位中文名和代码（从穴位统计表查找）
        apt_row = df_acupoints[df_acupoints["英文名"] == acupoint]
        name_cn = apt_row["中文名"].values[0] if len(apt_row) and pd.notna(apt_row["中文名"].values[0]) else ""
        code = apt_row["代码"].values[0] if len(apt_row) and pd.notna(apt_row["代码"].values[0]) else ""

        apt_info = get_acupoint_info(acupoint)
        results.append({
            "acupoint": acupoint,
            "name_cn": name_cn,
            "code": code,
            "cooccur": cooccur,
            "jaccard": jaccard,
            "gene_total": gene_total,
            "apt_total": apt_total,
            "article_count": len(article_list),
            "articles": article_list,
            "meridian": apt_info.get("meridian", ""),
            "apt_function": apt_info.get("function", ""),
        })

    return results


def get_network_data(min_cooccur=2, acupoints_filter=None, genes_filter=None, disease_filter=None, limit=300):
    """
    获取筛选后的网络数据
    acupoints_filter: 穴位英文名列表
    genes_filter: 基因符号列表
    disease_filter: 疾病名称
    """
    # 基础数据：从预计算共现网络开始（包含所有共现次数>=2的对）
    net = df_network[df_network["共现文献数"] >= min_cooccur].copy()

    if acupoints_filter:
        net = net[net["穴位"].isin(acupoints_filter)]

    if genes_filter:
        net = net[net["基因"].isin(genes_filter)]

    if disease_filter:
        # 找出该疾病的所有PMID
        disease_pmids = disease_pmids_map.get(disease_filter, set())
        if not disease_pmids:
            return []
        # 筛选出与该疾病相关的共现对
        # 从关联表中找出这些PMID对应的穴位-基因组合
        aa_d = df_aa[df_aa["PMID"].isin(disease_pmids)]
        ag_d = df_ag[df_ag["PMID"].isin(disease_pmids)]
        # 合并得到穴位-基因对
        merged = aa_d.merge(ag_d, on="PMID", how="inner")
        merged = merged.groupby(["穴位英文名", "基因符号"]).size().reset_index(name="count")
        merged = merged[merged["count"] >= min_cooccur]
        merged.columns = ["穴位", "基因", "共现文献数"]
        # 再应用穴位/基因筛选
        if acupoints_filter:
            merged = merged[merged["穴位"].isin(acupoints_filter)]
        if genes_filter:
            merged = merged[merged["基因"].isin(genes_filter)]
        net = merged

    net = net.nlargest(limit, "共现文献数")
    return [{"acupoint": r["穴位"], "gene": r["基因"], "count": int(r["共现文献数"])} for _, r in net.iterrows()]


def get_articles(page=1, per_page=20, search=""):
    """文献列表（支持搜索和分页）"""
    df = df_articles.copy()
    if search:
        mask = (
            df["标题"].astype(str).str.contains(search, case=False, na=False) |
            df["摘要"].astype(str).str.contains(search, case=False, na=False) |
            df["作者"].astype(str).str.contains(search, case=False, na=False)
        )
        df = df[mask]

    total = len(df)
    start = (page - 1) * per_page
    end = start + per_page
    df_page = df.iloc[start:end]

    articles = []
    for _, r in df_page.iterrows():
        authors = r["作者"] if pd.notna(r["作者"]) else ""
        articles.append({
            "id": int(r["ID"]),
            "pmid": r["PMID"],
            "title": r["标题"],
            "journal": r["期刊"] if pd.notna(r["期刊"]) else "",
            "year": int(r["年份"]) if pd.notna(r["年份"]) else None,
            "authors": authors[:80] + "..." if len(str(authors)) > 80 else authors,
            "doi": r["DOI"] if pd.notna(r["DOI"]) else "",
            "has_fulltext": bool(r["有全文"]) if pd.notna(r["有全文"]) else False,
        })

    return {"articles": articles, "total": total, "page": page, "per_page": per_page}


def get_article_detail(article_id):
    """文献详情"""
    row = df_articles[df_articles["ID"] == article_id]
    if row.empty:
        return None

    r = row.iloc[0]
    pmid = r["PMID"]

    article = {
        "id": int(r["ID"]),
        "pmid": pmid,
        "pmcid": r["PMCID"] if pd.notna(r["PMCID"]) else "",
        "doi": r["DOI"] if pd.notna(r["DOI"]) else "",
        "title": r["标题"],
        "abstract": r["摘要"] if pd.notna(r["摘要"]) else "",
        "journal": r["期刊"] if pd.notna(r["期刊"]) else "",
        "year": int(r["年份"]) if pd.notna(r["年份"]) else None,
        "authors": r["作者"] if pd.notna(r["作者"]) else "",
        "keywords": r["关键词"] if pd.notna(r["关键词"]) else "",
        "mesh": r["MeSH词"] if pd.notna(r["MeSH词"]) else "",
        "language": r["语言"] if pd.notna(r["语言"]) else "",
        "article_type": r["文章类型"] if pd.notna(r["文章类型"]) else "",
        "has_fulltext": bool(r["有全文"]) if pd.notna(r["有全文"]) else False,
        "search_strategy": r["检索策略"] if pd.notna(r["检索策略"]) else "",
    }

    # 关联基因
    genes = df_ag[df_ag["PMID"] == pmid]["基因符号"].unique().tolist()
    article["genes"] = genes

    # 关联穴位
    apt_rows = df_aa[df_aa["PMID"] == pmid][["穴位英文名", "代码"]].drop_duplicates()
    article["acupoints"] = [{"name": r["穴位英文名"], "code": r["代码"] if pd.notna(r["代码"]) else ""} for _, r in apt_rows.iterrows()]

    # 关联疾病
    diseases = df_ad[df_ad["PMID"] == pmid]["疾病"].unique().tolist()
    article["diseases"] = diseases

    return article


def get_journal_stats(limit=20):
    """期刊统计"""
    counts = df_articles[df_articles["期刊"].notna() & (df_articles["期刊"] != "")].groupby("期刊").size().reset_index(name="count")
    counts = counts.nlargest(limit, "count")
    return [{"journal": r["期刊"], "count": int(r["count"])} for _, r in counts.iterrows()]


# ========================================================================
# HTTP Handler
# ========================================================================

class Handler(BaseHTTPRequestHandler):
    def do_GET(self):
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path
        query = urllib.parse.parse_qs(parsed.query)

        try:
            # 页面路由
            if path == "/" or path == "/index.html":
                self._render_index()
            elif path == "/acupoint-query":
                self._render_acupoint_query(query)
            elif path == "/gene-query":
                self._render_gene_query(query)
            elif path == "/network":
                self._render_network(query)
            elif path == "/articles":
                page = int(query.get("page", ["1"])[0])
                search = query.get("q", [""])[0]
                self._render_articles(page, search)
            elif path.startswith("/article/"):
                article_id = int(path.split("/")[-1])
                self._render_article_detail(article_id)
            elif path == "/about":
                self._render_about()

            # API路由
            elif path == "/api/stats":
                self._json_response(get_stats())
            elif path == "/api/year_data":
                self._json_response(get_year_data())
            elif path == "/api/top_genes":
                limit = int(query.get("limit", ["20"])[0])
                self._json_response(get_top_genes(limit))
            elif path == "/api/top_acupoints":
                limit = int(query.get("limit", ["20"])[0])
                self._json_response(get_top_acupoints(limit))
            elif path == "/api/top_diseases":
                limit = int(query.get("limit", ["15"])[0])
                self._json_response(get_top_diseases(limit))
            elif path == "/api/journal_stats":
                limit = int(query.get("limit", ["20"])[0])
                self._json_response(get_journal_stats(limit))
            elif path == "/api/acupoint_list":
                self._json_response(get_acupoint_list())
            elif path == "/api/gene_list":
                self._json_response(get_gene_list())
            elif path == "/api/disease_list":
                self._json_response(get_disease_list())
            elif path == "/api/acupoint_genes":
                apt = query.get("acupoint", [""])[0]
                self._json_response(query_acupoint_genes(apt))
            elif path == "/api/gene_acupoints":
                gene = query.get("gene", [""])[0]
                self._json_response(query_gene_acupoints(gene))
            elif path == "/api/network":
                min_cooccur = int(query.get("min_cooccur", ["2"])[0])
                acupoints_filter = query.get("acupoints", [""])[0].split(",") if query.get("acupoints", [""])[0] else None
                genes_filter = query.get("genes", [""])[0].split(",") if query.get("genes", [""])[0] else None
                disease_filter = query.get("disease", [""])[0] or None
                self._json_response(get_network_data(min_cooccur, acupoints_filter, genes_filter, disease_filter))
            else:
                self._send_404()
        except Exception as e:
            self._send_error(str(e))

    # ---------- 页面渲染 ----------

    def _render_index(self):
        template = jenv.get_template("index.html")
        html = template.render(stats=get_stats())
        self._send_html(html)

    def _render_acupoint_query(self, query):
        template = jenv.get_template("acupoint_query.html")
        acupoint = query.get("acupoint", [""])[0]
        results = query_acupoint_genes(acupoint) if acupoint else []
        acupoint_list = get_acupoint_list()
        html = template.render(
            acupoint=acupoint,
            results=results,
            acupoint_list=acupoint_list
        )
        self._send_html(html)

    def _render_gene_query(self, query):
        template = jenv.get_template("gene_query.html")
        gene = query.get("gene", [""])[0]
        results = query_gene_acupoints(gene) if gene else []
        gene_list = get_gene_list()
        html = template.render(
            gene=gene,
            results=results,
            gene_list=gene_list
        )
        self._send_html(html)

    def _render_network(self, query):
        template = jenv.get_template("network.html")
        acupoint_list = get_acupoint_list()
        disease_list = get_disease_list()
        html = template.render(
            acupoint_list=acupoint_list,
            disease_list=disease_list
        )
        self._send_html(html)

    def _render_articles(self, page, search):
        template = jenv.get_template("articles.html")
        data = get_articles(page=page, search=search)
        html = template.render(
            articles=data["articles"],
            total=data["total"],
            page=data["page"],
            per_page=data["per_page"],
            search=search,
        )
        self._send_html(html)

    def _render_article_detail(self, article_id):
        template = jenv.get_template("article_detail.html")
        article = get_article_detail(article_id)
        if not article:
            self._send_404()
            return
        html = template.render(article=article)
        self._send_html(html)

    def _render_about(self):
        template = jenv.get_template("about.html")
        html = template.render()
        self._send_html(html)

    # ---------- HTTP 工具 ----------

    def _send_html(self, html):
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.end_headers()
        self.wfile.write(html.encode("utf-8"))

    def _json_response(self, data):
        self.send_response(200)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.end_headers()
        self.wfile.write(json.dumps(data, ensure_ascii=False, default=str).encode("utf-8"))

    def _send_404(self):
        self.send_response(404)
        self.send_header("Content-Type", "text/plain")
        self.end_headers()
        self.wfile.write(b"Not Found")

    def _send_error(self, msg):
        self.send_response(500)
        self.send_header("Content-Type", "text/plain")
        self.end_headers()
        self.wfile.write(f"Error: {msg}".encode("utf-8"))

    def log_message(self, format, *args):
        pass


def run_server(port=8080):
    server = HTTPServer(("0.0.0.0", port), Handler)
    print(f"\n🚀 针灸穴位-基因关联数据库分析平台 v2.1 已启动!")
    print(f"   数据源: {XLSX_PATH}")
    print(f"   本地访问: http://localhost:{port}")
    print(f"   按 Ctrl+C 停止\n")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\n服务器已停止")
        server.shutdown()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=8080)
    args = parser.parse_args()
    run_server(args.port)
