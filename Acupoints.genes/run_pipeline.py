#!/usr/bin/env python3
"""
针灸穴位-基因关联数据库
主运行脚本：文献采集 Pipeline

使用方法:
    python run_pipeline.py --init        # 初始化数据库
    python run_pipeline.py --fetch       # 执行文献采集
    python run_pipeline.py --stats       # 查看统计
    python run_pipeline.py --all         # 初始化+采集+统计
"""
import argparse
import sys
import os
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

from database.models import init_db, ArticleDB
from fetcher.pubmed_fetcher import PubMedFetcher
from config import SEARCH_STRATEGIES, DB_PATH


def init_database():
    """初始化数据库"""
    print("=" * 60)
    print("步骤 1/3: 初始化数据库")
    print("=" * 60)
    init_db(DB_PATH)
    print()


def fetch_articles(target_total: int = 300):
    """执行文献采集"""
    print("=" * 60)
    print(f"步骤 2/3: 文献采集 (目标: {target_total} 篇)")
    print("=" * 60)

    db = ArticleDB(DB_PATH)
    fetcher = PubMedFetcher()

    current_count = db.get_article_count()
    print(f"[INFO] 数据库当前已有: {current_count} 篇")

    if current_count >= target_total:
        print(f"[INFO] 已达到目标数量 ({current_count} >= {target_total})，跳过采集")
        return

    remaining = target_total - current_count
    total_new = 0

    for strategy in SEARCH_STRATEGIES:
        current_count = db.get_article_count()
        if current_count >= target_total:
            break

        print(f"\n[STRATEGY] {strategy['name']}: {strategy['query'][:60]}...")
        per_strategy_limit = min(remaining // 3 + 20, 100)  # 每个策略最多100篇

        try:
            result = fetcher.fetch_all_for_strategy(
                query=strategy["query"],
                strategy_name=strategy["name"],
                db=db,
                max_results=per_strategy_limit,
                skip_pdf=True
            )
            total_new += result["new"]
            print(f"[RESULT] 本策略新增: {result['new']}, 累计: {current_count + result['new']}")
        except Exception as e:
            print(f"[ERROR] 策略执行失败: {e}")
            continue

        # 避免请求过快，策略间休息
        time.sleep(2)

    final_count = db.get_article_count()
    print(f"\n[SUMMARY] 采集完成！新增: {final_count - current_count + total_new - (final_count - current_count)}")
    print(f"[SUMMARY] 数据库总量: {final_count} 篇")
    print()


def show_stats():
    """显示数据库统计"""
    print("=" * 60)
    print("步骤 3/3: 数据库统计")
    print("=" * 60)

    db = ArticleDB(DB_PATH)
    stats = db.get_stats()

    print(f"\n📊 数据库统计:")
    print(f"  - 文章总数: {stats['articles']}")
    print(f"  - 基因条目: {stats['genes']}")
    print(f"  - 穴位条目: {stats['acupoints']}")
    print(f"  - 疾病条目: {stats['diseases']}")

    print(f"\n📅 年份分布:")
    for year, count in stats.get("year_distribution", [])[-15:]:
        bar = "█" * (count // 2 + 1)
        print(f"  {year}: {bar} {count}")

    print()


def main():
    parser = argparse.ArgumentParser(description="针灸穴位-基因关联数据库 Pipeline")
    parser.add_argument("--init", action="store_true", help="初始化数据库")
    parser.add_argument("--fetch", action="store_true", help="执行文献采集")
    parser.add_argument("--stats", action="store_true", help="查看统计")
    parser.add_argument("--all", action="store_true", help="执行全部步骤")
    parser.add_argument("--target", type=int, default=300, help="目标文献数量 (默认300)")
    args = parser.parse_args()

    if args.all:
        args.init = True
        args.fetch = True
        args.stats = True

    if not (args.init or args.fetch or args.stats):
        parser.print_help()
        return

    if args.init:
        init_database()

    if args.fetch:
        fetch_articles(target_total=args.target)

    if args.stats:
        show_stats()


if __name__ == "__main__":
    main()
