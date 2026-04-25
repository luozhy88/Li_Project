"""
PubMed文献检索与下载模块
使用NCBI E-utilities API
"""
import time
import json
import os
import re
import xml.etree.ElementTree as ET
from urllib import request, parse, error
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import (
    NCBI_BASE_URL, NCBI_EMAIL, NCBI_API_KEY, REQUEST_DELAY,
    BATCH_SIZE, XML_DIR, PAPERS_DIR
)


class PubMedFetcher:
    """PubMed文献检索与元数据下载器"""

    def __init__(self, email: str = None, api_key: str = None):
        self.email = email or NCBI_EMAIL
        self.api_key = api_key or NCBI_API_KEY
        self.delay = REQUEST_DELAY
        self.last_request_time = 0
        os.makedirs(XML_DIR, exist_ok=True)
        os.makedirs(PAPERS_DIR, exist_ok=True)

    def _rate_limit(self):
        """请求限速"""
        elapsed = time.time() - self.last_request_time
        if elapsed < self.delay:
            time.sleep(self.delay - elapsed)
        self.last_request_time = time.time()

    def _get(self, url: str, retries: int = 3) -> str:
        """带重试的GET请求"""
        self._rate_limit()
        for attempt in range(retries):
            try:
                req = request.Request(url, headers={"User-Agent": f"AcupointGeneDB/1.0 ({self.email})"})
                with request.urlopen(req, timeout=15) as response:
                    return response.read().decode("utf-8")
            except error.HTTPError as e:
                if e.code == 429 or e.code >= 500:
                    wait = 2 ** attempt
                    print(f"  [WARN] HTTP {e.code}, 等待 {wait}s 后重试...")
                    time.sleep(wait)
                else:
                    raise
            except Exception as e:
                wait = 2 ** attempt
                print(f"  [WARN] 请求异常: {e}, 等待 {wait}s 后重试...")
                time.sleep(wait)
        raise Exception(f"请求失败: {url}")

    def search(self, query: str, retmax: int = 10000) -> list:
        """
        使用esearch检索PMID列表
        返回: PMID列表
        """
        print(f"[SEARCH] 检索: {query[:80]}...")
        params = {
            "db": "pubmed",
            "term": query,
            "retmode": "json",
            "retmax": retmax,
            "email": self.email,
        }
        if self.api_key:
            params["api_key"] = self.api_key

        url = f"{NCBI_BASE_URL}/esearch.fcgi?" + parse.urlencode(params)
        data = self._get(url)
        result = json.loads(data)
        pmids = result.get("esearchresult", {}).get("idlist", [])
        total = int(result.get("esearchresult", {}).get("count", 0))
        print(f"[SEARCH] 找到 {total} 篇, 本次获取 {len(pmids)} 篇")
        return pmids, total

    def fetch_details(self, pmids: list) -> list:
        """
        使用efetch获取文献详情 (XML -> 解析为dict)
        返回: 文章dict列表
        """
        if not pmids:
            return []
        print(f"[FETCH] 下载 {len(pmids)} 篇文献详情...")
        params = {
            "db": "pubmed",
            "id": ",".join(pmids),
            "retmode": "xml",
            "email": self.email,
        }
        if self.api_key:
            params["api_key"] = self.api_key

        url = f"{NCBI_BASE_URL}/efetch.fcgi?" + parse.urlencode(params)
        xml_data = self._get(url)

        # 保存原始XML
        timestamp = int(time.time())
        xml_path = os.path.join(XML_DIR, f"batch_{timestamp}.xml")
        with open(xml_path, "w", encoding="utf-8") as f:
            f.write(xml_data)

        articles = self._parse_pubmed_xml(xml_data)
        print(f"[FETCH] 成功解析 {len(articles)} 篇")
        return articles

    def _parse_pubmed_xml(self, xml_data: str) -> list:
        """解析PubMed XML为结构化数据"""
        articles = []
        try:
            root = ET.fromstring(xml_data.encode("utf-8"))
        except ET.ParseError as e:
            print(f"[ERROR] XML解析失败: {e}")
            return articles

        for article_elem in root.findall(".//PubmedArticle"):
            article = {}
            # PMID
            pmid_elem = article_elem.find(".//PMID")
            article["pmid"] = pmid_elem.text if pmid_elem is not None else None

            # DOI
            doi_elem = article_elem.find(".//ArticleId[@IdType='doi']")
            article["doi"] = doi_elem.text if doi_elem is not None else None

            # PMCID
            pmcid_elem = article_elem.find(".//ArticleId[@IdType='pmc']")
            article["pmcid"] = pmcid_elem.text if pmcid_elem is not None else None

            # 标题
            title_elem = article_elem.find(".//ArticleTitle")
            article["title"] = "".join(title_elem.itertext()) if title_elem is not None else ""

            # 摘要
            abstract_elem = article_elem.find(".//Abstract")
            if abstract_elem is not None:
                abstract_parts = []
                for abs_text in abstract_elem.findall("AbstractText"):
                    label = abs_text.get("Label", "")
                    text = "".join(abs_text.itertext())
                    if label:
                        abstract_parts.append(f"{label}: {text}")
                    else:
                        abstract_parts.append(text)
                article["abstract"] = "\n".join(abstract_parts)
            else:
                article["abstract"] = ""

            # 期刊
            journal_elem = article_elem.find(".//Journal/Title")
            article["journal"] = journal_elem.text if journal_elem is not None else ""

            # 日期
            pubdate = ""
            year_elem = article_elem.find(".//JournalIssue/PubDate/Year")
            if year_elem is not None:
                pubdate = year_elem.text
                article["year"] = int(year_elem.text) if year_elem.text.isdigit() else None
            else:
                medlinedate = article_elem.find(".//JournalIssue/PubDate/MedlineDate")
                if medlinedate is not None:
                    pubdate = medlinedate.text
                    # 提取年份
                    m = re.search(r'(\d{4})', pubdate)
                    article["year"] = int(m.group(1)) if m else None
                else:
                    article["year"] = None
            article["pubdate"] = pubdate

            # 作者
            authors = []
            affils = []
            author_list = article_elem.find(".//AuthorList")
            if author_list is not None:
                for author in author_list.findall("Author"):
                    last = author.find("LastName")
                    fore = author.find("ForeName")
                    if last is not None:
                        name = f"{last.text} {fore.text}" if fore is not None else last.text
                        authors.append(name)
                    #  affiliations
                    for affil in author.findall(".//Affiliation"):
                        if affil.text:
                            affils.append(affil.text)
            article["authors"] = "; ".join(authors)
            article["affiliations"] = "; ".join(list(set(affils)))

            # MeSH terms
            mesh_terms = []
            mesh_list = article_elem.find(".//MeshHeadingList")
            if mesh_list is not None:
                for mesh in mesh_list.findall("MeshHeading"):
                    desc = mesh.find("DescriptorName")
                    if desc is not None:
                        mesh_terms.append(desc.text)
            article["mesh_terms"] = "; ".join(mesh_terms)

            # Keywords
            keywords = []
            kw_list = article_elem.find(".//KeywordList")
            if kw_list is not None:
                for kw in kw_list.findall("Keyword"):
                    keywords.append("".join(kw.itertext()))
            article["keywords"] = "; ".join(keywords)

            # 语言
            lang_elem = article_elem.find(".//Language")
            article["language"] = lang_elem.text if lang_elem is not None else ""

            # 文章类型
            types = []
            for pt in article_elem.findall(".//PublicationType"):
                types.append(pt.text)
            article["article_type"] = "; ".join(types)

            # 是否有全文（PMCID存在则大概率有）
            article["has_fulltext"] = 1 if article.get("pmcid") else 0

            # XML路径
            article["xml_path"] = ""
            article["pdf_path"] = ""

            articles.append(article)

        return articles

    def fetch_pmc_pdf(self, pmcid: str, save_dir: str = None) -> str:
        """
        尝试下载PMC全文PDF（仅开放获取文章）
        返回: 保存路径或None
        """
        if not pmcid:
            return None
        save_dir = save_dir or PAPERS_DIR
        save_path = os.path.join(save_dir, f"{pmcid}.pdf")
        if os.path.exists(save_path) and os.path.getsize(save_path) > 1000:
            return save_path

        # PMC PDF链接格式
        url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/pdf/"
        try:
            self._rate_limit()
            req = request.Request(url, headers={"User-Agent": f"AcupointGeneDB/1.0 ({self.email})"})
            with request.urlopen(req, timeout=15) as response:
                data = response.read()
                if len(data) > 1000:
                    with open(save_path, "wb") as f:
                        f.write(data)
                    print(f"  [PDF] 下载成功: {pmcid}.pdf")
                    return save_path
        except Exception as e:
            print(f"  [PDF] 下载失败 {pmcid}: {e}")
        return None

    def fetch_all_for_strategy(self, query: str, strategy_name: str, db, max_results: int = None, skip_pdf: bool = False) -> dict:
        """
        执行完整的检索->下载->入库流程
        返回统计信息
        """
        pmids, total = self.search(query)
        if not pmids:
            db.log_fetch(strategy_name, query, total, 0, "no_results", "未找到文献")
            return {"total": 0, "new": 0, "duplicate": 0, "pmids": []}

        # 去重：检查数据库中已有的PMID
        existing_pmids = set(db.get_all_pmids())
        new_pmids = [p for p in pmids if p not in existing_pmids]
        duplicate_count = len(pmids) - len(new_pmids)

        print(f"[FETCH] 新文献: {len(new_pmids)}, 已存在: {duplicate_count}")

        if not new_pmids:
            db.log_fetch(strategy_name, query, total, 0, "all_duplicate", f"全部{duplicate_count}篇已存在")
            return {"total": total, "new": 0, "duplicate": duplicate_count, "pmids": pmids}

        # 限制数量
        if max_results and len(new_pmids) > max_results:
            new_pmids = new_pmids[:max_results]

        # 分批下载
        fetched_count = 0
        for i in range(0, len(new_pmids), BATCH_SIZE):
            batch = new_pmids[i:i+BATCH_SIZE]
            try:
                articles = self.fetch_details(batch)
                for art in articles:
                    art["search_strategy"] = strategy_name
                    db.insert_article(art)
                    fetched_count += 1
                    # 尝试下载PDF
                    if not skip_pdf and art.get("pmcid"):
                        self.fetch_pmc_pdf(art["pmcid"])
                print(f"  [PROGRESS] {fetched_count}/{len(new_pmids)}")
            except Exception as e:
                print(f"  [ERROR] 批次下载失败: {e}")

        db.log_fetch(strategy_name, query, total, fetched_count, "success",
                     f"duplicate={duplicate_count}")
        return {
            "total": total,
            "new": fetched_count,
            "duplicate": duplicate_count,
            "pmids": pmids
        }


if __name__ == "__main__":
    fetcher = PubMedFetcher()
    # 简单测试
    test_pmids, test_total = fetcher.search("acupuncture gene expression")
    print(f"测试检索: {test_total} 篇")
