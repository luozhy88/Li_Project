import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = ['Noto Sans CJK JP', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 数据源1: 图实际使用的
df1 = pd.read_csv('../SRP151288_HC_vs_TC.Genus/H20_feature_importance_heatmap_all.csv')
df1 = df1[df1['variable'].str.contains('Lactococcus')][['model','Importance']].set_index('model')

# 数据源2: 当前目录的CSV
df2 = pd.read_csv('H20_feature_importance_data.csv')
df2 = df2[df2['OTU'].str.contains('Lactococcus')][['DRF','GBM','GLM','DeepLearning']].T
df2.index.name = 'model'
df2.columns = ['Importance']

models = ['DRF', 'GBM', 'GLM', 'DeepLearning']
colors = ['#E74C3C', '#3498DB', '#2ECC71', '#9B59B6']

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# 左图: 数据源1
ax1 = axes[0]
bars1 = ax1.bar(models, [df1.loc[m, 'Importance'] for m in models], color=colors, edgecolor='black')
ax1.set_ylim(0, 1.0)
ax1.set_ylabel('Importance (scaled)', fontsize=12)
ax1.set_title('数据源1: H20_feature_importance_heatmap_all.csv\n(图 circular_bar_Genus_top50_grouped_from_csv 实际使用)', fontsize=10)
for bar, m in zip(bars1, models):
    h = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2, h + 0.02, f'{h:.3f}', ha='center', va='bottom', fontsize=11, weight='bold')

# 右图: 数据源2
ax2 = axes[1]
bars2 = ax2.bar(models, [df2.loc[m, 'Importance'] for m in models], color=colors, edgecolor='black')
ax2.set_ylim(0, 1.0)
ax2.set_ylabel('Importance', fontsize=12)
ax2.set_title('数据源2: H20_feature_importance_data.csv\n(当前目录下的文件)', fontsize=10)
for bar, m in zip(bars2, models):
    h = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2, h + 0.02, f'{h:.3f}', ha='center', va='bottom', fontsize=11, weight='bold')

plt.suptitle('Lactococcus 在两个数据源中的对比', fontsize=14, weight='bold')
plt.tight_layout()
fig.savefig('Lactococcus_data_comparison.pdf', dpi=300, bbox_inches='tight', facecolor='white')
fig.savefig('Lactococcus_data_comparison.svg', dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print('Saved: Lactococcus_data_comparison.pdf / .svg')
