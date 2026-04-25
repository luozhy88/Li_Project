"""
批量PDF下载工具：根据数据库中的PMCID下载开放获取全文
"""
import os
import sys
import sqlite3
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from fetcher.pubmed_fetcher import PubMedFetcher
from config import DB_PATH, PAPERS_DIR


def batch_download_pdfs(db_path: str = None, limit: int = None, delay: float = 1.0):
    """
    批量下载数据库中有PMCID但还没有PDF的文章
    
    Args:
        db_path: 数据库路径
        limit: 最多下载数量
        delay: 请求间隔（秒）
    """
    db_path = db_path or DB_PATH
    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    # 查询有PMCID但没有PDF的文章
    c.execute(
        """SELECT id, pmid, pmcid FROM articles 
           WHERE pmcid IS NOT NULL AND pmcid != '' 
           AND (pdf_path IS NULL OR pdf_path = '' OR NOT EXISTS (
               SELECT 1 FROM (SELECT 1) WHERE ? || '/' || pmcid || '.pdf' = pdf_path
           ))
           ORDER BY year DESC""",
        (PAPERS_DIR,)
    )
    rows = c.fetchall()
    conn.close()

    if limit:
        rows = rows[:limit]

    print(f"[PDF] 找到 {len(rows)} 篇待下载PDF")

    fetcher = PubMedFetcher()
    success = 0
    failed = 0

    for i, (article_id, pmid, pmcid) in enumerate(rows, 1):
        print(f"  [{i}/{len(rows)}] PMID:{pmid} PMCID:{pmcid} ...", end=" ")
        try:
            path = fetcher.fetch_pmc_pdf(pmcid)
            if path:
                # 更新数据库
                conn = sqlite3.connect(db_path)
                conn.execute("UPDATE articles SET pdf_path=?, has_fulltext=1 WHERE id=?", (path, article_id))
                conn.commit()
                conn.close()
                success += 1
                print("OK")
            else:
                failed += 1
                print("NOT_AVAILABLE")
        except Exception as e:
            failed += 1
            print(f"ERROR: {e}")

        time.sleep(delay)

    print(f"\n[PDF] 下载完成: 成功 {success}, 失败 {failed}")


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="批量下载PMC PDF")
    parser.add_argument("--limit", type=int, default=None, help="最多下载数量")
    parser.add_argument("--delay", type=float, default=1.0, help="请求间隔(秒)")
    args = parser.parse_args()
    batch_download_pdfs(limit=args.limit, delay=args.delay)
