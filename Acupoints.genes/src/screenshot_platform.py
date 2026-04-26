#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
截图分析平台关键页面 (使用Selenium + 系统Chromium)
"""

import os
import time
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By

OUTPUT_DIR = 'output/paper_figures'
os.makedirs(OUTPUT_DIR, exist_ok=True)

BASE_URL = 'http://ag.db.qyjc.top'


def get_driver():
    chrome_options = Options()
    chrome_options.add_argument('--headless')
    chrome_options.add_argument('--no-sandbox')
    chrome_options.add_argument('--disable-dev-shm-usage')
    chrome_options.add_argument('--window-size=1920,1080')
    # 使用系统chromium
    chrome_options.binary_location = '/usr/bin/chromium-browser'
    driver = webdriver.Chrome(options=chrome_options)
    return driver


def screenshot_pages():
    driver = get_driver()

    try:
        # 1. 首页 Dashboard
        print("正在截图首页 Dashboard...")
        driver.get(f'{BASE_URL}/')
        time.sleep(3)
        driver.save_screenshot(f'{OUTPUT_DIR}/fig1_dashboard.png')
        print(f"  已保存: {OUTPUT_DIR}/fig1_dashboard.png")

        # 2. 共现网络页面
        print("正在截图共现网络...")
        driver.get(f'{BASE_URL}/network')
        time.sleep(4)
        driver.save_screenshot(f'{OUTPUT_DIR}/fig2_network.png')
        print(f"  已保存: {OUTPUT_DIR}/fig2_network.png")

        # 3. 穴位查基因页面 - 查询足三里
        print("正在截图穴位查基因 (足三里)...")
        driver.get(f'{BASE_URL}/acupoint-query')
        time.sleep(2)
        try:
            input_box = driver.find_element(By.CSS_SELECTOR, 'input[type="text"]')
            input_box.clear()
            input_box.send_keys('足三里')
            submit_btn = driver.find_element(By.CSS_SELECTOR, 'button[type="submit"]')
            submit_btn.click()
            time.sleep(3)
        except Exception as e:
            print(f"  自动查询失败: {e}, 截图当前状态")
        driver.save_screenshot(f'{OUTPUT_DIR}/fig3_acupoint_query.png')
        print(f"  已保存: {OUTPUT_DIR}/fig3_acupoint_query.png")

        # 4. 基因查穴位页面 - 查询TNF
        print("正在截图基因查穴位 (TNF)...")
        driver.get(f'{BASE_URL}/gene-query')
        time.sleep(2)
        try:
            input_box = driver.find_element(By.CSS_SELECTOR, 'input[type="text"]')
            input_box.clear()
            input_box.send_keys('TNF')
            submit_btn = driver.find_element(By.CSS_SELECTOR, 'button[type="submit"]')
            submit_btn.click()
            time.sleep(3)
        except Exception as e:
            print(f"  自动查询失败: {e}, 截图当前状态")
        driver.save_screenshot(f'{OUTPUT_DIR}/fig4_gene_query.png')
        print(f"  已保存: {OUTPUT_DIR}/fig4_gene_query.png")

        # 5. 文献列表页面
        print("正在截图文献浏览页面...")
        driver.get(f'{BASE_URL}/articles')
        time.sleep(3)
        driver.save_screenshot(f'{OUTPUT_DIR}/fig5_articles.png')
        print(f"  已保存: {OUTPUT_DIR}/fig5_articles.png")

    finally:
        driver.quit()

    print("\n所有截图完成！")


if __name__ == '__main__':
    screenshot_pages()
