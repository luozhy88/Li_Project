#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Figure 5: transcriptome validation composite (volcano + enrichment forest plot)."""
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path

sup = Path('edit.Paper/output/20261024_25_transcriptome_validation')
out = Path('edit.Paper/output/20261024_03_support/20261024_Figure5.png')

de = pd.read_csv(sup / '20261024_03_all_measured_gene_DE.csv')
cov = pd.read_csv(sup / '20261024_05_candidate_coverage.csv')
enr = pd.read_csv(sup / '20261024_04_candidate_set_enrichment.csv')
ci = pd.read_csv(sup / '20261025_figure5_fisher_ci.csv')

sets = {r['set']: set(r['measured_genes'].split(';')) for _, r in cov.iterrows()}
colors = {'ZBM': '#E69F00', 'XKC': '#009E73', 'XYT': '#7A5195'}  # orange / teal / purple, match Figures 2-3

# point estimates must match the archived enrichment table
for _, r in ci.iterrows():
    ref = enr[(enr['set'] == r['set']) & (enr['background'] == r['bg'])].iloc[0]
    assert abs(r['or'] - ref['odds_ratio']) < 1e-6, (r['set'], r['bg'])
print('OR point estimates match archived enrichment table')

fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.2), gridspec_kw={'width_ratios': [1.15, 1]})

# --- Panel A: volcano ---
ax = axes[0]
x = de['logFC'].to_numpy(); y = -np.log10(de['adj.P.Val'].to_numpy())
genes = de['gene'].to_numpy()
ax.scatter(x, y, s=2, c='#BBBBBB', alpha=0.35, linewidths=0, rasterized=True)
for s in ['XYT', 'XKC', 'ZBM']:  # smallest set drawn last (on top)
    m = np.isin(genes, list(sets[s]))
    ax.scatter(x[m], y[m], s=4, c=colors[s], alpha=0.75, linewidths=0,
               label=f'{s} (n={len(sets[s])})', rasterized=True)
ax.axvline(1, ls='--', lw=0.8, c='k'); ax.axvline(-1, ls='--', lw=0.8, c='k')
ax.axhline(-np.log10(0.05), ls='--', lw=0.8, c='k')
ax.set_xlabel('log2 fold change (primary PTC vs paired normal thyroid)')
ax.set_ylabel('−log10(FDR)')
ax.legend(loc='upper left', fontsize=8, frameon=False, markerscale=2.5)
ax.set_title('A', loc='left', fontsize=13, fontweight='bold')

# --- Panel B: forest plot ---
ax = axes[1]
rows = []
for s in ['XYT', 'XKC', 'ZBM']:
    for bg_name, bg_lab in [('all_measured', 'all measured genes'), ('disease_measured', 'disease-measured genes')]:
        r = ci[(ci['set'] == s) & (ci['bg'] == bg_name)].iloc[0]
        fdr = enr[(enr['set'] == s) & (enr['background'] == bg_name)].iloc[0]['FDR_6tests']
        rows.append((s, bg_lab, r['or'], r['lo'], r['hi'], fdr))
ypos = np.arange(len(rows))[::-1]
for (s, bg_lab, orv, lo, hi, fdr), yi in zip(rows, ypos):
    ax.plot([lo, hi], [yi, yi], color=colors[s], lw=1.6, solid_capstyle='round')
    ax.plot(orv, yi, 'o', color=colors[s], ms=6)
    ax.text(6.0, yi, f'{orv:.2f} ({lo:.2f}–{hi:.2f})  FDR={fdr:.1e}'.replace('e-0', 'e−').replace('e-', 'e−'),
            va='center', fontsize=7.5)
ax.axvline(1, ls='--', lw=0.8, c='k')
ax.set_xscale('log'); ax.set_xlim(0.25, 5.2)
ax.set_xlabel('Odds ratio (95% CI), log scale')
ax.set_yticks(ypos)
ax.set_yticklabels([f'{s} vs {bg}' for s, bg, *_ in rows], fontsize=8.5)
ax.set_title('B', loc='left', fontsize=13, fontweight='bold')
for a in axes:
    a.spines[['top', 'right']].set_visible(False)

fig.tight_layout()
fig.savefig(out, dpi=300)
print('saved', out)
