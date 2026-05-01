# ThyGlycoPortal - Cloudflare Worker Edition

Thyroid Cancer Glycosylation Interactive Analysis Platform deployed on Cloudflare Worker.

## Features

| Route | Description |
|-------|-------------|
| `/` | Overview Dashboard: stats cards + year trend + cancer/method distribution + biomarker AUC |
| `/literature` | Literature table with PMID links |
| `/biomarkers` | Biomarker performance table + AUC comparison chart |
| `/glycans` | Glycan browser with cross-study comparison |
| `/enzymes` | Enzyme-glycan regulatory network |
| `/tcga` | TCGA glycogene expression (Bones et al. 2018) |
| `/diagnostic` | Interactive BN score & recurrence ratio calculators |
| `/nomogram` | PTMC LNM risk nomogram (demonstration) |
| `/api/*` | RESTful JSON API |

## Tech Stack

- **Runtime**: Cloudflare Worker (V8 Isolate)
- **Data Source**: `output/thyroid_glyco_db.sqlite` -> inline JSON (~25KB)
- **Charts**: ECharts 5 (CDN)
- **Style**: Tailwind CSS (CDN)
- **Frontend**: Pure HTML + inline JS (no framework)

## Deploy

```bash
# 1. Login to Cloudflare
npx wrangler login

# 2. Deploy
cd cloudflare-test-worker
npx wrangler deploy

# 3. Access URL will be printed
```

## Local Dev

```bash
cd cloudflare-test-worker
npm install
npm run dev
# Open http://localhost:8787
```

## Rebuild from Database

When the SQLite database is updated:

```bash
cd /path/to/Tangji
python3 cloudflare-test-worker/rebuild_index.py
npx wrangler deploy
```

## Project Structure

```
cloudflare-test-worker/
├── src/
│   └── index.js          # Worker main file (inline data + pages + API)
├── rebuild_index.py      # Rebuild src/index.js from SQLite
├── wrangler.toml         # Wrangler config
├── package.json          # Node.js config
└── README.md             # This file
```

## Data Sources

- Zhang ZJ et al. Front Oncol. 2021 (PMID: 34476207) — IgG BN diagnostic biomarker
- Kudelka MR et al. Cancer Med. 2023 (PMID: 36437732) — Serum recurrence prediction
- PTMC Nomogram. Front Oncol. 2022 (PMC9497917) — Lymph node metastasis risk
- Bones J et al. Cancers. 2018 (PMC11727208) — TCGA RNAseq glycogene expression
- Wu et al. Front Endocrinol. 2021 (PMC8267918) — Plasma glycomics BTN/TC

All data are real published summary statistics. No simulated data.
