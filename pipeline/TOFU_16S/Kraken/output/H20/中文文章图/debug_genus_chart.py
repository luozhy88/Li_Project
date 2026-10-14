import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = ['Noto Sans CJK JP', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

DPI = 600
models = ['DRF', 'GBM', 'GLM', 'DeepLearning']
colors_models = {
    'DRF': '#E74C3C',
    'GBM': '#3498DB',
    'GLM': '#2ECC71',
    'DeepLearning': '#9B59B6'
}

df = pd.read_csv('../SRP151288_HC_vs_TC.Genus/H20_feature_importance_heatmap_all.csv')
df_wide = df.pivot(index='variable', columns='model', values='Importance')
df_wide = df_wide.reset_index()
df_wide.columns.name = None
for m in models:
    if m not in df_wide.columns:
        df_wide[m] = np.nan
df_wide = df_wide[['variable'] + models]
df_wide = df_wide.dropna(subset=models, how='all').reset_index(drop=True)
df_wide['MeanImportance'] = df_wide[models].mean(axis=1)
df_wide['Label'] = df_wide['variable'].str.replace(r'^OTU\d+_g__', '', regex=True)

df_top = df_wide.sort_values('MeanImportance', ascending=True).tail(50).reset_index(drop=True)
n_bars = len(df_top)
n_models = len(models)

print(f"Total bars: {n_bars}")
print(f"Lactococcus index in df_top: {df_top[df_top['Label'] == 'Lactococcus'].index[0]}")
print(df_top[df_top['Label'] == 'Lactococcus'][['Label'] + models].to_string(index=False))
print()

fig = plt.figure(figsize=(20, 20))
ax = fig.add_subplot(111, polar=True)

total_slots = n_bars * n_models
theta_all = np.linspace(0, 2 * np.pi, total_slots, endpoint=False)
width = 2 * np.pi / total_slots * 0.82

bars_per_model = []
for m_idx, model in enumerate(models):
    indices = np.arange(m_idx, total_slots, n_models)
    theta_model = theta_all[indices]
    radii = df_top[model].fillna(0).values
    bars = ax.bar(theta_model, radii, width=width, bottom=0.0,
           color=colors_models[model], edgecolor='white', linewidth=0.25,
           alpha=0.85, label=model)
    bars_per_model.append(bars)

# 检查 Lactococcus 对应的 bars
lacto_idx = df_top[df_top['Label'] == 'Lactococcus'].index[0]
print(f"\nLactococcus bar index: {lacto_idx}")
for m_idx, model in enumerate(models):
    bar = bars_per_model[m_idx][lacto_idx]
    print(f"{model}: height={bar.get_height():.6f}, theta={bar.get_x():.6f}")

plt.close()
