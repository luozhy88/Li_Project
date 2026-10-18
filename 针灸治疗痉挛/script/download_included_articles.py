from __future__ import annotations

import csv
import re
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urljoin
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[1]
CSV_PATH = ROOT / "output/02_PubMed初步分析/01_逐篇筛选与书目.csv"
OUT = ROOT / "output/05_纳入文章"
PUBMED_OUT = OUT / "01_PubMed"
WOS_OUT = OUT / "02_Web_of_Science_待提供导出"


def clean_title(title: str) -> str:
    # Keep the article title recognizable; only replace filesystem-forbidden characters.
    title = re.sub(r"[\\/:*?\"<>|\x00-\x1f]", "／", title).strip().rstrip(".")
    return title or "untitled"


def fetch(url: str) -> tuple[bytes, str] | None:
    req = Request(url, headers={"User-Agent": "Mozilla/5.0 (research article retrieval)"})
    try:
        with urlopen(req, timeout=45) as response:
            return response.read(), response.headers.get_content_type()
    except (HTTPError, URLError, TimeoutError):
        return None


def write_local_abstract(row: dict[str, str], path: Path) -> None:
    abstract = row.get("abstract", "").strip() or "摘要未收录于当前 PubMed 记录。"
    doi = row.get("doi", "").strip()
    doi_line = f'<p><b>DOI：</b><a href="https://doi.org/{doi}">{doi}</a></p>' if doi else ""
    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{row['title']}</title>
<style>body{{font-family:Arial,sans-serif;max-width:900px;margin:40px auto;padding:0 24px;line-height:1.65;color:#222}}h1{{line-height:1.25}}.meta{{color:#555;border-bottom:1px solid #ddd;padding-bottom:1em}}.abstract{{white-space:pre-wrap}}</style></head>
<body><h1>{row['title']}</h1><div class="meta"><p><b>PMID：</b><a href="{row['url']}">{row['pmid']}</a></p>{doi_line}</div>
<h2>Abstract</h2><div class="abstract">{abstract}</div>
</body></html>"""
    path.write_text(html, encoding="utf-8")


def save_article(row: dict[str, str]) -> dict[str, str]:
    title = clean_title(row["title"])
    pmid = row["pmid"]
    pmc = row["pmc"].strip()
    doi = row["doi"].strip()
    candidates: list[str] = []
    if pmc:
        candidates.extend([
            f"https://pmc.ncbi.nlm.nih.gov/articles/{pmc}/pdf/",
            f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc}/pdf/",
        ])
    if doi:
        candidates.append(f"https://doi.org/{doi}")
        if doi.startswith("10.3389/"):
            journal = doi.removeprefix("10.3389/").split(".", 1)[0]
            journal_slug = {"fneur": "neurology"}.get(journal, journal)
            candidates.append(f"https://www.frontiersin.org/journals/{journal_slug}/articles/{doi}/pdf")

    for url in candidates:
        result = fetch(url)
        if result and result[1] == "application/pdf" and result[0][:4] == b"%PDF":
            path = PUBMED_OUT / f"{title}.pdf"
            path.write_bytes(result[0])
            return {"pmid": pmid, "title": row["title"], "status": "full_text_pdf", "file": path.name, "source": url}

    # PMC sometimes serves the complete article as HTML while withholding the PDF.
    if pmc:
        pmc_url = f"https://pmc.ncbi.nlm.nih.gov/articles/{pmc}/"
        result = fetch(pmc_url)
        if result and result[1] == "text/html" and len(result[0]) > 10000:
            html = result[0].decode("utf-8", errors="ignore")
            if "recaptcha" in html.lower() or "challengepage" in html.lower():
                result = None
            else:
                pdf_links = re.findall(r"href=[\"']([^\"']+\.pdf(?:\?[^\"']*)?)[\"']", html, re.I)
                for link in pdf_links:
                    pdf_result = fetch(urljoin(pmc_url, link))
                    if pdf_result and pdf_result[1] == "application/pdf" and pdf_result[0][:4] == b"%PDF":
                        path = PUBMED_OUT / f"{title}.pdf"
                        path.write_bytes(pdf_result[0])
                        return {"pmid": pmid, "title": row["title"], "status": "full_text_pdf", "file": path.name, "source": urljoin(pmc_url, link)}
                path = PUBMED_OUT / f"{title}.html"
                path.write_bytes(result[0])
                return {"pmid": pmid, "title": row["title"], "status": "full_text_html", "file": path.name, "source": pmc_url}

    # Full text was not publicly downloadable; preserve the PubMed abstract page.
    result = fetch(row["url"])
    path = PUBMED_OUT / f"{title}.html"
    # PubMed may return a small anti-bot/challenge document that appears blank locally.
    if result and len(result[0]) > 10000:
        pubmed_html = result[0].decode("utf-8", errors="ignore").lower()
        if "recaptcha" not in pubmed_html and "challengepage" not in pubmed_html:
            path.write_bytes(result[0])
            return {"pmid": pmid, "title": row["title"], "status": "abstract_page", "file": path.name, "source": row["url"]}
    write_local_abstract(row, path)
    return {"pmid": pmid, "title": row["title"], "status": "local_abstract", "file": path.name, "source": row["url"]}


def main() -> None:
    PUBMED_OUT.mkdir(parents=True, exist_ok=True)
    WOS_OUT.mkdir(parents=True, exist_ok=True)
    with CSV_PATH.open(encoding="utf-8-sig", newline="") as f:
        rows = [row for row in csv.DictReader(f) if row["relevance_decision"] == "include"]

    results = []
    for row in rows:
        results.append(save_article(row))
        time.sleep(0.4)

    with (OUT / "下载清单.csv").open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["pmid", "title", "status", "file", "source"])
        writer.writeheader()
        writer.writerows(results)

    (OUT / "README.md").write_text(
        "# 纳入文章\n\n"
        "本目录按数据库分组保存纳入文章。当前项目中仅发现 PubMed 检索结果；\n"
        "Web of Science 的原始导出文件尚未提供，因此 `02_Web_of_Science_待提供导出` 暂为空，\n"
        "未将 PubMed 记录冒充为 WoS 记录。\n\n"
        "PubMed 子目录中的 PDF/全文 HTML 为自动获取到的公开全文；无法取得公开全文或遇到验证码的记录，保存为可离线打开的本地摘要 HTML。\n"
        "具体下载状态见 `下载清单.csv`。\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
