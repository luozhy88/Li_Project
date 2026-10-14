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

# 读取原始数据
df_raw = pd.read_csv('../SRP151288_HC_vs_TC.Genus/H20_feature_importance_heatmap_all.csv')

# === 第一步：按 scaled Importance 排序取 Top 20（与圆形图一致） ===
df_scaled = df_raw.pivot(index='variable', columns='model', values='Importance').reset_index()
df_scaled.columns.name = None
df_scaled['MeanImportance'] = df_scaled[models].mean(axis=1)
df_scaled['Label'] = df_scaled['variable'].str.replace(r'^OTU\d+_g__', '', regex=True)
top20_vars = df_scaled.sort_values('MeanImportance', ascending=False).head(20)['variable'].tolist()

# === 第二步：对这20个菌属，提取 percentage 值 ===
df_pct = df_raw[df_raw['variable'].isin(top20_vars)].pivot(
    index='variable', columns='model', values='percentage').reset_index()
df_pct.columns.name = None

# 合并标签
df_pct = df_pct.merge(df_scaled[['variable', 'Label']], on='variable', how='left')

# 按 scaled MeanImportance 降序排列，然后反转使Top1在最上面
df_pct = df_pct.set_index('variable').reindex(top20_vars).reset_index()
df_pct = df_pct.iloc[::-1].reset_index(drop=True)

# 截断特别长的名字
for i, label in enumerate(df_pct['Label']):
    if len(label) > 16:
        df_pct.at[i, 'Label'] = label[:7] + '...' + label[-4:]

fig, ax = plt.subplots(figsize=(12, 14))

n_bars = len(df_pct)
n_models = len(models)
bar_height = 0.14
group_gap = 0.18

for m_idx, model in enumerate(models):
    y_positions = []
    for i in range(n_bars):
        group_center = i * (n_models * bar_height + group_gap)
        y = group_center + (n_models - 1 - m_idx) * bar_height
        y_positions.append(y)
    
    values = df_pct[model].values
    ax.barh(y_positions, values, height=bar_height,
            color=colors_models[model], edgecolor='white', linewidth=0.3,
            alpha=0.9, label=model_labels[model])

# y轴标签
y_tick_positions = []
for i in range(n_bars):
    group_center = i * (n_models * bar_height + group_gap)
    y_tick_positions.append(group_center + (n_models - 1) * bar_height / 2 + bar_height / 2)

ax.set_yticks(y_tick_positions)
ax.set_yticklabels(df_pct['Label'].values, fontsize=13)
# 给y轴顶端和底端留一点边距，确保最顶端/最底端刻度线不被裁剪
ax.margins(y=0.03)

# 横轴：percentage（各模型内总和为1）
max_val = df_pct[models].values.max()
ax.set_xlabel('特征重要性', fontsize=14)
ax.set_xlim(0, max_val * 1.2)

# 刻度朝里，只在底部和左侧显示刻度线（默认行为）
ax.tick_params(axis='both', direction='in', length=6, width=1)

# 竖向网格线
ax.xaxis.grid(True, linestyle='--', alpha=0.4)
ax.set_axisbelow(True)

# 只显示底部和左侧边框
for spine in ['top', 'right']:
    ax.spines[spine].set_visible(False)

# 图例：右下角，汉字
ax.legend(loc='lower right', fontsize=12, frameon=False)

# 标题
ax.set_title('Top20菌属特征重要性', fontsize=16, pad=15, weight='bold')

plt.tight_layout()
fig.savefig('barplot_top20_genus_percentage.pdf', dpi=DPI, bbox_inches='tight', facecolor='white')
fig.savefig('barplot_top20_genus_percentage.svg', dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close()

print('Saved: barplot_top20_genus_percentage.pdf / .svg')
print('Top 20 菌属（与圆形图一致）:', list(df_pct['Label'].values))
