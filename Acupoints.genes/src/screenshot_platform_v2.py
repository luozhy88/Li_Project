#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
截图分析平台关键页面 v2
"""

import os
import time
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait, Select
from selenium.webdriver.support import expected_conditions as EC

OUTPUT_DIR = 'output/platform_screenshots'
os.makedirs(OUTPUT_DIR, exist_ok=True)

BASE_URL = 'http://ag.db.qyjc.top'


def get_driver():
    chrome_options = Options()
    chrome_options.add_argument('--headless')
    chrome_options.add_argument('--no-sandbox')
    chrome_options.add_argument('--disable-dev-shm-usage')
    chrome_options.add_argument('--window-size=1920,1080')
    chrome_options.add_argument('--disable-gpu')
    chrome_options.binary_location = '/usr/bin/chromium-browser'
    driver = webdriver.Chrome(options=chrome_options)
    return driver


def screenshot_pages():
    driver = get_driver()
    wait = WebDriverWait(driver, 10)

    try:
        # 1. 首页 Dashboard
        print("正在截图首页 Dashboard...")
        driver.get(f'{BASE_URL}/')
        wait.until(EC.presence_of_element_located((By.CSS_SELECTOR, '.stat-card')))
        time.sleep(2)
        driver.save_screenshot(f'{OUTPUT_DIR}/fig1_dashboard.png')
        print(f"  已保存: fig1_dashboard.png")

        # 2. 共现网络页面
        print("正在截图共现网络...")
        driver.get(f'{BASE_URL}/network')
        wait.until(EC.presence_of_element_located((By.ID, 'network-chart')))
        time.sleep(3)
        driver.save_screenshot(f'{OUTPUT_DIR}/fig2_network.png')
        print(f"  已保存: fig2_network.png")

        # 3. 穴位查基因页面 - 查询 Zusanli
        print("正在截图穴位查基因 (Zusanli)...")
        driver.get(f'{BASE_URL}/acupoint-query')
        wait.until(EC.presence_of_element_located((By.CSS_SELECTOR, 'select[name="acupoint"]')))
        time.sleep(1)
        try:
            select = Select(driver.find_element(By.CSS_SELECTOR, 'select[name="acupoint"]'))
            select.select_by_value('Zusanli')
            submit_btn = driver.find_element(By.CSS_SELECTOR, 'button[type="submit"]')
            submit_btn.click()
            wait.until(EC.presence_of_element_located((By.ID, 'gene-bar-chart')))
            time.sleep(2)
        except Exception as e:
            print(f"  自动查询失败: {e}, 截图当前状态")
        driver.save_screenshot(f'{OUTPUT_DIR}/fig3_acupoint_query.png')
        print(f"  已保存: fig3_acupoint_query.png")

        # 4. 基因查穴位页面 - 查询 TNF
        print("正在截图基因查穴位 (TNF)...")
        driver.get(f'{BASE_URL}/gene-query')
        wait.until(EC.presence_of_element_located((By.CSS_SELECTOR, 'input[name="gene"]')))
        time.sleep(1)
        try:
            input_box = driver.find_element(By.CSS_SELECTOR, 'input[name="gene"]')
            input_box.clear()
            input_box.send_keys('TNF')
            submit_btn = driver.find_element(By.CSS_SELECTOR, 'button[type="submit"]')
            submit_btn.click()
            wait.until(EC.presence_of_element_located((By.ID, 'apt-bar-chart')))
            time.sleep(2)
        except Exception as e:
            print(f"  自动查询失败: {e}, 截图当前状态")
        driver.save_screenshot(f'{OUTPUT_DIR}/fig4_gene_query.png')
        print(f"  已保存: fig4_gene_query.png")

        # 5. 文献列表页面
        print("正在截图文献浏览页面...")
        driver.get(f'{BASE_URL}/articles')
        wait.until(EC.presence_of_element_located((By.CSS_SELECTOR, '.data-table')))
        time.sleep(2)
        driver.save_screenshot(f'{OUTPUT_DIR}/fig5_articles.png')
        print(f"  已保存: fig5_articles.png")

        # 6. 算法说明页面
        print("正在截图算法说明页面...")
        driver.get(f'{BASE_URL}/about')
        wait.until(EC.presence_of_element_located((By.CSS_SELECTOR, '.card')))
        time.sleep(2)
        driver.save_screenshot(f'{OUTPUT_DIR}/fig6_about.png')
        print(f"  已保存: fig6_about.png")

    except Exception as e:
        print(f"截图过程中出错: {e}")
        import traceback
        traceback.print_exc()
    finally:
        driver.quit()

    print("\n所有截图完成！")


if __name__ == '__main__':
    screenshot_pages()
