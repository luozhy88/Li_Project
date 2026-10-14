import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

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
    'DRF': '分布式随机森林',
    'GBM': '梯度提升机',
    'GLM': '广义线性模型',
    'DeepLearning': '深度学习'
}


def load_and_pivot(csv_path):
    """读取长格式CSV，转换为宽格式（OTU为行，模型为列）"""
    df = pd.read_csv(csv_path)
    # 确保model列的值正确
    df_wide = df.pivot(index='variable', columns='model', values='Importance')
    df_wide = df_wide.reset_index()
    df_wide.columns.name = None
    # 确保4个模型列都存在
    for m in models:
        if m not in df_wide.columns:
            df_wide[m] = np.nan
    df_wide = df_wide[['variable'] + models]
    df_wide = df_wide.dropna(subset=models, how='all').reset_index(drop=True)
    df_wide['MeanImportance'] = df_wide[models].mean(axis=1)
    df_wide['MaxImportance'] = df_wide[models].max(axis=1)
    # 提取标签名
    if '_g__' in df_wide['variable'].iloc[0]:
        df_wide['Label'] = df_wide['variable'].str.replace(r'^OTU\d+_g__', '', regex=True)
    elif '_s__' in df_wide['variable'].iloc[0]:
        df_wide['Label'] = df_wide['variable'].str.replace(r'^OTU\d+_s__', '', regex=True)
    else:
        df_wide['Label'] = df_wide['variable']
    return df_wide


def draw_top50_grouped(df_in, title, filename_base, show_every=1):
    """绘制Top 50分组环形图"""
    df_top = df_in.sort_values('MeanImportance', ascending=True).tail(50).reset_index(drop=True)
    n_bars = len(df_top)
    n_models = len(models)

    fig = plt.figure(figsize=(20, 20))
    ax = fig.add_subplot(111, polar=True)

    total_slots = n_bars * n_models
    theta_all = np.linspace(0, 2 * np.pi, total_slots, endpoint=False)
    width = 2 * np.pi / total_slots * 0.82

    for m_idx, model in enumerate(models):
        indices = np.arange(m_idx, total_slots, n_models)
        theta_model = theta_all[indices]
        radii = df_top[model].fillna(0).values
        ax.bar(theta_model, radii, width=width, bottom=0.0,
               color=colors_models[model], edgecolor='white', linewidth=0.25,
               alpha=0.85, label=model)

    # 标签：每隔一个显示，减少密度
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

    ax.set_xticks([])
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=10, color='gray')
    ax.set_ylim(0, 1.05)
    ax.spines['polar'].set_visible(False)

    legend_elements = [Patch(facecolor=colors_models[m], edgecolor='white', label=model_labels[m]) for m in models]
    ax.legend(handles=legend_elements, loc='upper right',
              bbox_to_anchor=(1.03, 1.01), fontsize=10, frameon=False)

    plt.title('', fontsize=16, pad=12, weight='bold')

    # 保存SVG和PDF
    for ext in ['.svg', '.pdf']:
        fname = filename_base + ext
        fig.savefig(fname, dpi=DPI, bbox_inches='tight', facecolor='white', pad_inches=0.08)
        print(f"Saved: {fname}")
    plt.close()


if __name__ == '__main__':
    # === 属级 ===
    df_genus = load_and_pivot('../SRP151288_HC_vs_TC.Genus/H20_feature_importance_heatmap_all.csv')
    print(f"Genus OTUs: {len(df_genus)}")
    draw_top50_grouped(
        df_genus,
        'Top 50 Genus OTUs Feature Importance\n(Grouped by 4 Models)',
        'circular_bar_Genus_top50_grouped_from_csv',
        show_every=1
    )

    # === 种级 ===
    df_species = load_and_pivot('../SRP151288_HC_vs_TC.Species/H20_feature_importance_heatmap_all.csv')
    print(f"Species OTUs: {len(df_species)}")
    draw_top50_grouped(
        df_species,
        'Top 50 Species OTUs Feature Importance\n(Grouped by 4 Models)',
        'circular_bar_species_top50_grouped_from_csv',
        show_every=1
    )

    print("\n全部完成！")
