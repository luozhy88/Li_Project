import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = ['Noto Sans CJK JP', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

DPI = 300
models = ['DRF', 'GBM', 'GLM', 'DeepLearning']
colors_models = {
    'DRF': '#E74C3C',
    'GBM': '#3498DB',
    'GLM': '#2ECC71',
    'DeepLearning': '#9B59B6'
}

model_labels = {
    'DRF': '分布式随机森林',
    'GBM': '梯度提升机',
    'GLM': '广义线性模型',
    'DeepLearning': '深度学习'
}

# 读取数据
df = pd.read_csv('../SRP151288_HC_vs_TC.Genus/H20_feature_importance_heatmap_all.csv')
df_wide = df.pivot(index='variable', columns='model', values='Importance').reset_index()
df_wide.columns.name = None
df_wide['MeanImportance'] = df_wide[models].mean(axis=1)
df_wide['Label'] = df_wide['variable'].str.replace(r'^OTU\d+_g__', '', regex=True)

# Top 20 降序后反转，Top1在最上面
df_top = df_wide.sort_values('MeanImportance', ascending=False).head(20).iloc[::-1].reset_index(drop=True)

# 截断特别长的名字
for i, label in enumerate(df_top['Label']):
    if len(label) > 16:
        df_top.at[i, 'Label'] = label[:7] + '...' + label[-4:]

fig, ax = plt.subplots(figsize=(14, 18))

n_bars = len(df_top)
n_models = len(models)
bar_height = 0.15
group_gap = 0.35

for m_idx, model in enumerate(models):
    y_positions = []
    for i in range(n_bars):
        group_center = i * (n_models * bar_height + group_gap)
        y = group_center + (n_models - 1 - m_idx) * bar_height
        y_positions.append(y)
    
    values = df_top[model].values
    ax.barh(y_positions, values, height=bar_height,
            color=colors_models[model], edgecolor='white', linewidth=0.3,
            alpha=0.9, label=model_labels[model])

# y轴标签
y_tick_positions = []
for i in range(n_bars):
    group_center = i * (n_models * bar_height + group_gap)
    y_tick_positions.append(group_center + (n_models - 1) * bar_height / 2 + bar_height / 2)

ax.set_yticks(y_tick_positions)
ax.set_yticklabels(df_top['Label'].values, fontsize=13)

# 横轴
ax.set_xlabel('标准化重要性', fontsize=14)
ax.set_xlim(0, 1.05)
ax.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])

# 刻度朝里
ax.tick_params(axis='both', direction='in', length=6, width=1)

# 竖向网格线
ax.xaxis.grid(True, linestyle='--', alpha=0.4)
ax.set_axisbelow(True)

# 只显示底部和左侧边框，隐藏顶部和右侧
for spine in ['top', 'right']:
    ax.spines[spine].set_visible(False)

# 图例：右上角，汉字
ax.legend(loc='upper right', fontsize=12, frameon=False)

# 标题
ax.set_title('Top20菌属特征重要性', fontsize=16, pad=15, weight='bold')

plt.tight_layout()
fig.savefig('barplot_top20_genus.pdf', dpi=DPI, bbox_inches='tight', facecolor='white')
fig.savefig('barplot_top20_genus.svg', dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close()

print('Saved: barplot_top20_genus.pdf / .svg')
