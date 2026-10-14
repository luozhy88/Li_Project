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
model_labels = {
    'DRF': 'DRF',
    'GBM': 'GBM',
    'GLM': 'GLM',
    'DeepLearning': 'DL'
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

# 标注所有标签
show_every = 2
label_indices = np.arange(1, total_slots, n_models * show_every)
radii_full = np.zeros(total_slots)
labels_full = np.empty(total_slots, dtype=object)
for i in range(n_bars):
    for j in range(n_models):
        slot = i * n_models + j
        radii_full[slot] = df_top[models].iloc[i].max()
        labels_full[slot] = df_top['Label'].iloc[i]

for li in label_indices:
    angle = theta_all[li]
    radius = radii_full[li] + 0.03
    label = labels_full[li]
    if len(label) > 22:
        label = label[:19] + '...'
    angle_deg = np.degrees(angle)
    if 90 < angle_deg <= 270:
        ha = 'right'
        text_angle = angle_deg + 180
    else:
        ha = 'left'
        text_angle = angle_deg
    ax.text(angle, radius, label,
            ha=ha, va='center', fontsize=7,
            rotation=text_angle - 90,
            rotation_mode='anchor',
            color='#111111', weight='bold')

# 在 Lactococcus 柱子上方标注数值
lacto_idx = df_top[df_top['Label'] == 'Lactococcus'].index[0]
print(f"Lactococcus index: {lacto_idx}")
for m_idx, model in enumerate(models):
    bar = bars_per_model[m_idx][lacto_idx]
    h = bar.get_height()
    x = bar.get_x() + bar.get_width()/2
    ax.text(x, h + 0.05, f'{model_labels[model]}:{h:.2f}', 
            ha='center', va='bottom', fontsize=6, color=colors_models[model],
            rotation=np.degrees(x)-90, rotation_mode='anchor')
    print(f"  {model}: {h:.6f}")

ax.set_xticks([])
ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax.set_yticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=10, color='gray')
ax.set_ylim(0, 1.05)
ax.spines['polar'].set_visible(False)
ax.legend(loc='upper right', bbox_to_anchor=(1.03, 1.01), fontsize=10, frameon=False)

fig.savefig('circular_bar_Genus_top50_Lactococcus_debug.svg', dpi=DPI, bbox_inches='tight', facecolor='white', pad_inches=0.08)
fig.savefig('circular_bar_Genus_top50_Lactococcus_debug.pdf', dpi=DPI, bbox_inches='tight', facecolor='white', pad_inches=0.08)
plt.close()
print("Saved debug figure.")
