"""
命令行分析工具：数据库统计、可视化、导出
"""
import sqlite3
import os
import sys
import json
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config import DB_PATH, OUTPUT_DIR


class CLIAnalyzer:
    """数据库分析器"""

    def __init__(self, db_path: str = None):
        self.db_path = db_path or DB_PATH
        self.conn = sqlite3.connect(self.db_path)
        self.conn.row_factory = sqlite3.Row

    def close(self):
        self.conn.close()

    def _query(self, sql, params=()):
        c = self.conn.cursor()
        c.execute(sql, params)
        return c.fetchall()

    def overview(self):
        """数据库概览"""
        print("\n" + "=" * 60)
        print("📊 数据库概览")
        print("=" * 60)

        counts = {
            "articles": self._query("SELECT COUNT(*) FROM articles")[0][0],
            "genes": self._query("SELECT COUNT(*) FROM genes")[0][0],
            "acupoints": self._query("SELECT COUNT(*) FROM acupoints")[0][0],
            "diseases": self._query("SELECT COUNT(*) FROM diseases")[0][0],
            "article_gene": self._query("SELECT COUNT(*) FROM article_gene")[0][0],
            "article_acupoint": self._query("SELECT COUNT(*) FROM article_acupoint")[0][0],
            "article_disease": self._query("SELECT COUNT(*) FROM article_disease")[0][0],
        }

        for k, v in counts.items():
            print(f"  {k:20s}: {v}")

        # 有全文的数量
        ft = self._query("SELECT COUNT(*) FROM articles WHERE has_fulltext=1")[0][0]
        print(f"  {'has_fulltext':20s}: {ft}")

        return counts

    def year_trend(self):
        """年份趋势分析"""
        print("\n" + "=" * 60)
        print("📅 发表年份趋势 (近15年)")
        print("=" * 60)

        rows = self._query(
            "SELECT year, COUNT(*) as cnt FROM articles WHERE year IS NOT NULL GROUP BY year ORDER BY year DESC LIMIT 15"
        )
        for row in rows:
            bar = "█" * min(row[1], 50)
            print(f"  {row[0]}: {bar} {row[1]}")
        return rows

    def top_journals(self, n: int = 15):
        """ Top期刊 """
        print("\n" + "=" * 60)
        print("📚 Top 期刊")
        print("=" * 60)

        rows = self._query(
            "SELECT journal, COUNT(*) as cnt FROM articles WHERE journal != '' GROUP BY journal ORDER BY cnt DESC LIMIT ?",
            (n,)
        )
        for i, row in enumerate(rows, 1):
            print(f"  {i:2d}. {row[0][:45]:45s} ({row[1]})")
        return rows

    def top_genes(self, n: int = 20):
        """高频基因"""
        print("\n" + "=" * 60)
        print("🧬 Top 高频基因")
        print("=" * 60)

        rows = self._query(
            """SELECT g.gene_symbol, COUNT(ag.article_id) as cnt
               FROM genes g JOIN article_gene ag ON g.id = ag.gene_id
               GROUP BY g.id ORDER BY cnt DESC LIMIT ?""",
            (n,)
        )
        for i, row in enumerate(rows, 1):
            print(f"  {i:2d}. {row[0]:15s} ({row[1]} 篇)")
        return rows

    def top_acupoints(self, n: int = 20):
        """高频穴位"""
        print("\n" + "=" * 60)
        print("📍 Top 高频穴位")
        print("=" * 60)

        rows = self._query(
            """SELECT a.name_en, COUNT(aa.article_id) as cnt
               FROM acupoints a JOIN article_acupoint aa ON a.id = aa.acupoint_id
               GROUP BY a.id ORDER BY cnt DESC LIMIT ?""",
            (n,)
        )
        for i, row in enumerate(rows, 1):
            print(f"  {i:2d}. {row[0]:15s} ({row[1]} 篇)")
        return rows

    def top_diseases(self, n: int = 15):
        """高频疾病"""
        print("\n" + "=" * 60)
        print("🏥 Top 高频疾病/表型")
        print("=" * 60)

        rows = self._query(
            """SELECT d.name, COUNT(ad.article_id) as cnt
               FROM diseases d JOIN article_disease ad ON d.id = ad.disease_id
               GROUP BY d.id ORDER BY cnt DESC LIMIT ?""",
            (n,)
        )
        for i, row in enumerate(rows, 1):
            print(f"  {i:2d}. {row[0]:30s} ({row[1]} 篇)")
        return rows

    def acupoint_gene_network(self, min_cooccur: int = 2):
        """穴位-基因共现网络"""
        print("\n" + "=" * 60)
        print("🔗 穴位-基因共现网络 (共现≥{}次)".format(min_cooccur))
        print("=" * 60)

        rows = self._query(
            """SELECT a.name_en, g.gene_symbol, COUNT(*) as cnt
               FROM article_acupoint aa
               JOIN article_gene ag ON aa.article_id = ag.article_id
               JOIN acupoints a ON aa.acupoint_id = a.id
               JOIN genes g ON ag.gene_id = g.id
               GROUP BY a.id, g.id HAVING cnt >= ?
               ORDER BY cnt DESC LIMIT 30""",
            (min_cooccur,)
        )
        for row in rows:
            print(f"  {row[0]:12s} ↔ {row[1]:12s} ({row[2]} 篇)")
        return rows

    def export_summary_json(self, filepath: str = None):
        """导出汇总JSON"""
        filepath = filepath or os.path.join(OUTPUT_DIR, "summary.json")
        os.makedirs(os.path.dirname(filepath), exist_ok=True)

        data = {
            "overview": {
                "articles": self._query("SELECT COUNT(*) FROM articles")[0][0],
                "genes": self._query("SELECT COUNT(*) FROM genes")[0][0],
                "acupoints": self._query("SELECT COUNT(*) FROM acupoints")[0][0],
                "diseases": self._query("SELECT COUNT(*) FROM diseases")[0][0],
            },
            "year_distribution": [
                dict(r) for r in self._query(
                    "SELECT year, COUNT(*) as count FROM articles WHERE year IS NOT NULL GROUP BY year ORDER BY year"
                )
            ],
            "top_genes": [
                {"gene": r[0], "count": r[1]} for r in self._query(
                    """SELECT g.gene_symbol, COUNT(*) as cnt FROM genes g
                       JOIN article_gene ag ON g.id = ag.gene_id
                       GROUP BY g.id ORDER BY cnt DESC LIMIT 50"""
                )
            ],
            "top_acupoints": [
                {"acupoint": r[0], "count": r[1]} for r in self._query(
                    """SELECT a.name_en, COUNT(*) as cnt FROM acupoints a
                       JOIN article_acupoint aa ON a.id = aa.acupoint_id
                       GROUP BY a.id ORDER BY cnt DESC LIMIT 50"""
                )
            ],
            "top_diseases": [
                {"disease": r[0], "count": r[1]} for r in self._query(
                    """SELECT d.name, COUNT(*) as cnt FROM diseases d
                       JOIN article_disease ad ON d.id = ad.disease_id
                       GROUP BY d.id ORDER BY cnt DESC LIMIT 30"""
                )
            ],
        }

        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False, indent=2)
        print(f"\n[EXPORT] 汇总JSON已导出: {filepath}")
        return data

    def export_csv(self, outdir: str = None):
        """导出CSV表格"""
        import csv
        outdir = outdir or OUTPUT_DIR
        os.makedirs(outdir, exist_ok=True)

        # 文章列表
        with open(os.path.join(outdir, "articles.csv"), "w", newline="", encoding="utf-8-sig") as f:
            writer = csv.writer(f)
            writer.writerow(["pmid", "title", "journal", "year", "authors", "doi", "has_fulltext"])
            for row in self._query("SELECT pmid, title, journal, year, authors, doi, has_fulltext FROM articles"):
                writer.writerow(row)
        print(f"[EXPORT] articles.csv -> {outdir}")

        # 基因列表
        with open(os.path.join(outdir, "genes.csv"), "w", newline="", encoding="utf-8-sig") as f:
            writer = csv.writer(f)
            writer.writerow(["gene_symbol", "article_count"])
            for row in self._query(
                """SELECT g.gene_symbol, COUNT(*) FROM genes g
                   JOIN article_gene ag ON g.id = ag.gene_id
                   GROUP BY g.id ORDER BY COUNT(*) DESC"""
            ):
                writer.writerow(row)
        print(f"[EXPORT] genes.csv -> {outdir}")

        # 穴位列表
        with open(os.path.join(outdir, "acupoints.csv"), "w", newline="", encoding="utf-8-sig") as f:
            writer = csv.writer(f)
            writer.writerow(["name_en", "article_count"])
            for row in self._query(
                """SELECT a.name_en, COUNT(*) FROM acupoints a
                   JOIN article_acupoint aa ON a.id = aa.acupoint_id
                   GROUP BY a.id ORDER BY COUNT(*) DESC"""
            ):
                writer.writerow(row)
        print(f"[EXPORT] acupoints.csv -> {outdir}")

    def full_report(self):
        """生成完整分析报告"""
        self.overview()
        self.year_trend()
        self.top_journals()
        self.top_acupoints()
        self.top_genes()
        self.top_diseases()
        self.acupoint_gene_network(min_cooccur=2)
        self.export_summary_json()
        self.export_csv()


def main():
    import argparse
    parser = argparse.ArgumentParser(description="针灸穴位-基因数据库分析工具")
    parser.add_argument("--overview", action="store_true", help="数据库概览")
    parser.add_argument("--year", action="store_true", help="年份趋势")
    parser.add_argument("--journals", action="store_true", help="Top期刊")
    parser.add_argument("--genes", action="store_true", help="Top基因")
    parser.add_argument("--acupoints", action="store_true", help="Top穴位")
    parser.add_argument("--diseases", action="store_true", help="Top疾病")
    parser.add_argument("--network", action="store_true", help="共现网络")
    parser.add_argument("--export", action="store_true", help="导出JSON和CSV")
    parser.add_argument("--all", action="store_true", help="执行全部分析")
    args = parser.parse_args()

    if not any([args.overview, args.year, args.journals, args.genes, args.acupoints,
                args.diseases, args.network, args.export, args.all]):
        args.all = True

    analyzer = CLIAnalyzer()

    if args.all:
        analyzer.full_report()
    else:
        if args.overview: analyzer.overview()
        if args.year: analyzer.year_trend()
        if args.journals: analyzer.top_journals()
        if args.genes: analyzer.top_genes()
        if args.acupoints: analyzer.top_acupoints()
        if args.diseases: analyzer.top_diseases()
        if args.network: analyzer.acupoint_gene_network()
        if args.export:
            analyzer.export_summary_json()
            analyzer.export_csv()

    analyzer.close()


if __name__ == "__main__":
    main()
