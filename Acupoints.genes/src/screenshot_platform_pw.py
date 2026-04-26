#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
截图分析平台关键页面 (Playwright + 系统 Chromium)
"""

import os
from playwright.sync_api import sync_playwright

OUTPUT_DIR = 'output/platform_screenshots'
os.makedirs(OUTPUT_DIR, exist_ok=True)

BASE_URL = 'http://ag.db.qyjc.top'


def screenshot_pages():
    with sync_playwright() as p:
        browser = p.chromium.launch(
            headless=True,
            executable_path='/usr/bin/chromium-browser',
            args=['--no-sandbox', '--disable-dev-shm-usage']
        )
        context = browser.new_context(viewport={'width': 1920, 'height': 1080})
        page = context.new_page()

        # 1. 首页 Dashboard
        print("正在截图首页 Dashboard...")
        page.goto(f'{BASE_URL}/', wait_until='domcontentloaded', timeout=60000)
        page.wait_for_timeout(3000)
        page.screenshot(path=f'{OUTPUT_DIR}/fig1_dashboard.png', full_page=False)
        print(f"  已保存: fig1_dashboard.png")

        # 2. 共现网络页面
        print("正在截图共现网络...")
        page.goto(f'{BASE_URL}/network', wait_until='domcontentloaded', timeout=60000)
        page.wait_for_timeout(4000)
        page.screenshot(path=f'{OUTPUT_DIR}/fig2_network.png', full_page=False)
        print(f"  已保存: fig2_network.png")

        # 3. 穴位查基因页面 - 查询 Zusanli
        print("正在截图穴位查基因 (Zusanli)...")
        page.goto(f'{BASE_URL}/acupoint-query', wait_until='domcontentloaded', timeout=60000)
        page.wait_for_timeout(1500)
        page.select_option('select[name="acupoint"]', 'Zusanli')
        page.click('button[type="submit"]')
        page.wait_for_selector('#gene-bar-chart', state='visible', timeout=15000)
        page.wait_for_timeout(2000)
        page.screenshot(path=f'{OUTPUT_DIR}/fig3_acupoint_query.png', full_page=False)
        print(f"  已保存: fig3_acupoint_query.png")

        # 4. 基因查穴位页面 - 查询 TNF
        print("正在截图基因查穴位 (TNF)...")
        page.goto(f'{BASE_URL}/gene-query', wait_until='domcontentloaded', timeout=60000)
        page.wait_for_timeout(1500)
        page.fill('input[name="gene"]', 'TNF')
        page.click('button[type="submit"]')
        page.wait_for_selector('#apt-bar-chart', state='visible', timeout=15000)
        page.wait_for_timeout(2000)
        page.screenshot(path=f'{OUTPUT_DIR}/fig4_gene_query.png', full_page=False)
        print(f"  已保存: fig4_gene_query.png")

        # 5. 文献列表页面
        print("正在截图文献浏览页面...")
        page.goto(f'{BASE_URL}/articles', wait_until='domcontentloaded', timeout=60000)
        page.wait_for_timeout(2000)
        page.screenshot(path=f'{OUTPUT_DIR}/fig5_articles.png', full_page=False)
        print(f"  已保存: fig5_articles.png")

        # 6. 算法说明页面
        print("正在截图算法说明页面...")
        page.goto(f'{BASE_URL}/about', wait_until='domcontentloaded', timeout=60000)
        page.wait_for_timeout(2000)
        page.screenshot(path=f'{OUTPUT_DIR}/fig6_about.png', full_page=False)
        print(f"  已保存: fig6_about.png")

        browser.close()

    print("\n所有截图完成！")


if __name__ == '__main__':
    screenshot_pages()
