# 针灸穴位-基因关联数据库 - Cloudflare Worker 版

基于 Cloudflare Worker 部署的全功能 Web 分析平台，数据全部内联在 Worker 脚本中，无需额外数据库或存储。

## 功能特性

| 页面 | 说明 |
|------|------|
| `/` | Dashboard：统计卡片 + 年份趋势图 + Top排行 |
| `/acupoint-query` | 穴位 → 基因查询（共现频率 + Jaccard + 生物功能 + 支持文献） |
| `/gene-query` | 基因 → 穴位查询（共现频率 + Jaccard + 归经与主治 + 支持文献） |
| `/network` | 共现网络可视化（力导向图 + 穴位/基因/疾病三维筛选） |
| `/articles` | 文献列表（分页 + 搜索） |
| `/article/:id` | 文献详情（含摘要、关联基因/穴位/疾病） |
| `/about` | 算法说明（共现频率、Jaccard公式、网络指标） |
| `/api/*` | RESTful API，返回 JSON 数据 |

## 技术栈

- **运行时**: Cloudflare Worker (V8 Isolate)
- **数据源**: `output/acupoint_gene_database.xlsx` → 内联 JSON (~470KB)
- **图表**: ECharts 5 (CDN)
- **前端**: 纯 HTML + 内联 CSS/JS（无框架依赖）

## 快速部署

### 方式一：Wrangler CLI（推荐）

```bash
# 1. 登录 Cloudflare
npx wrangler login

# 2. 进入项目目录
cd cloudflare-test-worker

# 3. 部署
npx wrangler deploy

# 4. 获取访问地址
# 输出类似：https://acupoint-gene-db.<your-subdomain>.workers.dev
```

### 方式二：Cloudflare Dashboard（手动）

1. 登录 [Cloudflare Dashboard](https://dash.cloudflare.com)
2. 进入 **Workers & Pages** → **创建 Worker**
3. 编辑代码，将 `src/index.js` 的全部内容粘贴进去
4. 点击 **部署**
5. 获得 `https://<name>.workers.dev` 访问地址

## 本地开发

```bash
cd cloudflare-test-worker
npm install
npm run dev
# 访问 http://localhost:8787
```

> 注：macOS 12.x 及以下版本可能遇到 wrangler 运行时兼容性问题，建议直接在 Dashboard 中部署。

## 数据更新

当 Excel 数据更新后，重新生成 Worker 代码：

```bash
cd /Users/apple/Desktop/github.local/Li_Project/Acupoints.genes
python3 src/web/export_worker_data.py   # 重新导出数据到 src/index.js
npx wrangler deploy
```

## 项目结构

```
cloudflare-test-worker/
├── src/
│   └── index.js          # Worker 主文件（含内联数据 + 查询逻辑 + HTML模板）
├── wrangler.toml         # Wrangler 配置
├── package.json          # 项目配置
└── README.md             # 本文件
```

## 限制说明

- **Worker 脚本大小**: ~520KB，在 Cloudflare Free Plan 的 1MB 限制内
- **文献摘要**: 为控制体积，摘要未内联存储；文献详情页提供 PubMed 外链
- **并发**: Cloudflare Worker 自动扩展，无并发限制
