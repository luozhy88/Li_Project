"""
数据库模型：针灸穴位-基因关联数据库 (SQLite)
"""
import sqlite3
from datetime import datetime
import os

# 数据库路径由外部传入或使用默认
DEFAULT_DB_PATH = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "output", "acupoint_gene.db")


def init_db(db_path: str = None):
    """初始化数据库，创建所有表"""
    db_path = db_path or DEFAULT_DB_PATH
    os.makedirs(os.path.dirname(db_path), exist_ok=True)
    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    # 文章主表
    c.execute('''
        CREATE TABLE IF NOT EXISTS articles (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            pmid TEXT UNIQUE,
            pmcid TEXT,
            doi TEXT,
            title TEXT,
            abstract TEXT,
            journal TEXT,
            pubdate TEXT,
            year INTEGER,
            authors TEXT,
            affiliations TEXT,
            keywords TEXT,
            mesh_terms TEXT,
            language TEXT,
            article_type TEXT,
            has_fulltext INTEGER DEFAULT 0,
            pdf_path TEXT,
            xml_path TEXT,
            search_strategy TEXT,
            created_at TEXT DEFAULT CURRENT_TIMESTAMP,
            updated_at TEXT DEFAULT CURRENT_TIMESTAMP
        )
    ''')

    # 基因表
    c.execute('''
        CREATE TABLE IF NOT EXISTS genes (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            gene_symbol TEXT,
            gene_name TEXT,
            entrez_id TEXT,
            uniprot_id TEXT,
            UNIQUE(gene_symbol)
        )
    ''')

    # 穴位表
    c.execute('''
        CREATE TABLE IF NOT EXISTS acupoints (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            name_en TEXT,
            name_cn TEXT,
            code TEXT,
            meridian TEXT,
            UNIQUE(code)
        )
    ''')

    # 文章-基因关联表
    c.execute('''
        CREATE TABLE IF NOT EXISTS article_gene (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            article_id INTEGER,
            gene_id INTEGER,
            regulation TEXT,  -- up/down/unknown
            evidence TEXT,    -- microarray/RNA-seq/qPCR/WB等
            context TEXT,     -- 组织/细胞类型
            FOREIGN KEY (article_id) REFERENCES articles(id),
            FOREIGN KEY (gene_id) REFERENCES genes(id),
            UNIQUE(article_id, gene_id)
        )
    ''')

    # 文章-穴位关联表
    c.execute('''
        CREATE TABLE IF NOT EXISTS article_acupoint (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            article_id INTEGER,
            acupoint_id INTEGER,
            stimulation TEXT,  -- 针刺/艾灸/电针/按压
            evidence TEXT,
            FOREIGN KEY (article_id) REFERENCES articles(id),
            FOREIGN KEY (acupoint_id) REFERENCES acupoints(id),
            UNIQUE(article_id, acupoint_id)
        )
    ''')

    # 疾病/表型表
    c.execute('''
        CREATE TABLE IF NOT EXISTS diseases (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            name TEXT,
            mesh_term TEXT,
            category TEXT,
            UNIQUE(name)
        )
    ''')

    # 文章-疾病关联表
    c.execute('''
        CREATE TABLE IF NOT EXISTS article_disease (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            article_id INTEGER,
            disease_id INTEGER,
            FOREIGN KEY (article_id) REFERENCES articles(id),
            FOREIGN KEY (disease_id) REFERENCES diseases(id),
            UNIQUE(article_id, disease_id)
        )
    ''')

    # 采集日志表
    c.execute('''
        CREATE TABLE IF NOT EXISTS fetch_log (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            strategy_name TEXT,
            query TEXT,
            total_count INTEGER,
            fetched_count INTEGER,
            status TEXT,
            message TEXT,
            created_at TEXT DEFAULT CURRENT_TIMESTAMP
        )
    ''')

    conn.commit()
    conn.close()
    print(f"[DB] 数据库初始化完成: {db_path}")


class ArticleDB:
    """数据库操作封装"""

    def __init__(self, db_path: str = None):
        self.db_path = db_path or DEFAULT_DB_PATH

    def _connect(self):
        return sqlite3.connect(self.db_path)

    def insert_article(self, article: dict) -> int:
        """插入文章，返回article_id"""
        conn = self._connect()
        c = conn.cursor()
        fields = [
            "pmid", "pmcid", "doi", "title", "abstract", "journal", "pubdate",
            "year", "authors", "affiliations", "keywords", "mesh_terms",
            "language", "article_type", "has_fulltext", "pdf_path", "xml_path",
            "search_strategy"
        ]
        placeholders = ",".join(["?"] * len(fields))
        sql = f"INSERT OR IGNORE INTO articles ({','.join(fields)}) VALUES ({placeholders})"
        values = [article.get(f) for f in fields]
        c.execute(sql, values)
        conn.commit()
        # 获取id
        c.execute("SELECT id FROM articles WHERE pmid=?", (article.get("pmid"),))
        row = c.fetchone()
        conn.close()
        return row[0] if row else None

    def get_article_count(self) -> int:
        conn = self._connect()
        c = conn.cursor()
        c.execute("SELECT COUNT(*) FROM articles")
        count = c.fetchone()[0]
        conn.close()
        return count

    def get_all_pmids(self) -> list:
        conn = self._connect()
        c = conn.cursor()
        c.execute("SELECT pmid FROM articles WHERE pmid IS NOT NULL")
        rows = c.fetchall()
        conn.close()
        return [r[0] for r in rows]

    def log_fetch(self, strategy_name: str, query: str, total_count: int, fetched_count: int, status: str, message: str = ""):
        conn = self._connect()
        c = conn.cursor()
        c.execute(
            "INSERT INTO fetch_log (strategy_name, query, total_count, fetched_count, status, message) VALUES (?,?,?,?,?,?)",
            (strategy_name, query, total_count, fetched_count, status, message)
        )
        conn.commit()
        conn.close()

    def get_stats(self) -> dict:
        """获取数据库统计信息"""
        conn = self._connect()
        c = conn.cursor()
        stats = {}
        c.execute("SELECT COUNT(*) FROM articles")
        stats["articles"] = c.fetchone()[0]
        c.execute("SELECT COUNT(*) FROM genes")
        stats["genes"] = c.fetchone()[0]
        c.execute("SELECT COUNT(*) FROM acupoints")
        stats["acupoints"] = c.fetchone()[0]
        c.execute("SELECT COUNT(*) FROM diseases")
        stats["diseases"] = c.fetchone()[0]
        c.execute("SELECT year, COUNT(*) FROM articles WHERE year IS NOT NULL GROUP BY year ORDER BY year")
        stats["year_distribution"] = c.fetchall()
        conn.close()
        return stats


if __name__ == "__main__":
    init_db()
