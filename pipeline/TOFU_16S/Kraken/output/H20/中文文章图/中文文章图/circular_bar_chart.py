import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os

# 设置中文字体 - 使用系统可用的 Noto Sans CJK JP
plt.rcParams['font.family'] = ['Noto Sans CJK JP', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 读取数据
df = pd.read_csv('H20_feature_importance_data.csv')
models = ['DRF', 'GBM', 'GLM', 'DeepLearning']

# 丢弃全为NaN的行
df = df.dropna(subset=models, how='all').reset_index(drop=True)

# 计算综合重要性（平均值和最大值）
df['MeanImportance'] = df[models].mean(axis=1)
df['MaxImportance'] = df[models].max(axis=1)

# 按平均重要性排序（升序，这样tail取的是最大的）
df_sorted = df.sort_values('MeanImportance', ascending=True).reset_index(drop=True)

# 提取OTU名称（去掉OTU999_前缀，只保留属名）
df_sorted['Label'] = df_sorted['OTU'].str.replace(r'^OTU\d+_', '', regex=True)

# 颜色方案
colors_models = {
    'DRF': '#E74C3C',
    'GBM': '#3498DB',
    'GLM': '#2ECC71',
    'DeepLearning': '#9B59B6'
}


def draw_circular_bar_single(df_in, value_col, title, filename, max_bars=195, label_top_n=40):
    """绘制单个环形条形图"""
    df_plot = df_in.tail(max_bars).reset_index(drop=True)
    n_bars = len(df_plot)

    fig = plt.figure(figsize=(14, 14))
    ax = fig.add_subplot(111, polar=True)

    # 角度
    theta = np.linspace(0, 2 * np.pi, n_bars, endpoint=False)
    width = 2 * np.pi / n_bars * 0.85

    # 半径
    radii = df_plot[value_col].values

    # 颜色：根据数值大小映射
    norm = plt.Normalize(radii.min(), radii.max())
    cmap = plt.cm.viridis
    bar_colors = cmap(norm(radii))

    # 绘制条形
    ax.bar(theta, radii, width=width, bottom=0.05, color=bar_colors, edgecolor='white', linewidth=0.3)

    # 只给最重要的前 label_top_n 个加标签
    label_indices = list(range(max(0, n_bars - label_top_n), n_bars))

    for idx in label_indices:
        angle = theta[idx]
        radius = radii[idx] + 0.08
        label = df_plot['Label'].iloc[idx]
        if len(label) > 20:
            label = label[:17] + '...'

        angle_deg = np.degrees(angle)
        if 90 < angle_deg <= 270:
            ha = 'right'
            angle_deg += 180
        else:
            ha = 'left'

        ax.text(angle, radius, label,
                ha=ha, va='center', fontsize=6,
                rotation=angle_deg - 90,
                rotation_mode='anchor',
                color='#333333')

    # 设置网格和刻度
    ax.set_xticks([])
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=9, color='gray')
    ax.set_ylim(0, 1.35)
    ax.spines['polar'].set_visible(False)

    # 添加颜色条（用inset_axes避免polar上的bug）
    cax = fig.add_axes([0.85, 0.15, 0.03, 0.2])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, cax=cax, label='Importance')

    plt.title(title, fontsize=16, pad=30, weight='bold')
    plt.savefig(filename, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {filename}")


def draw_circular_bar_facet(df_in, title_prefix, filename):
    """绘制2x2分面环形图，每个模型一个环"""
    fig, axes = plt.subplots(2, 2, figsize=(18, 18), subplot_kw=dict(polar=True))
    axes = axes.flatten()

    n_bars = len(df_in)
    theta = np.linspace(0, 2 * np.pi, n_bars, endpoint=False)
    width = 2 * np.pi / n_bars * 0.9

    for idx, model in enumerate(models):
        ax = axes[idx]
        radii = df_in[model].values

        norm = plt.Normalize(0, 1.0)
        cmap = plt.cm.plasma
        bar_colors = cmap(norm(radii))

        ax.bar(theta, radii, width=width, bottom=0.02, color=bar_colors, edgecolor='white', linewidth=0.2)

        # 只给Top 30加标签
        top_n = 30
        top_indices = np.argsort(radii)[-top_n:]

        for ti in top_indices:
            angle = theta[ti]
            radius = radii[ti] + 0.06
            label = df_in['Label'].iloc[ti]
            if len(label) > 18:
                label = label[:15] + '...'

            angle_deg = np.degrees(angle)
            if 90 < angle_deg <= 270:
                ha = 'right'
                angle_deg += 180
            else:
                ha = 'left'

            ax.text(angle, radius, label,
                    ha=ha, va='center', fontsize=4.5,
                    rotation=angle_deg - 90,
                    rotation_mode='anchor',
                    color='#222222')

        ax.set_xticks([])
        ax.set_yticks([0.25, 0.5, 0.75, 1.0])
        ax.set_yticklabels(['0.25', '0.5', '0.75', '1.0'], fontsize=7, color='gray')
        ax.set_ylim(0, 1.25)
        ax.spines['polar'].set_visible(False)
        ax.set_title(model, fontsize=14, pad=15, weight='bold', color=colors_models[model])

    plt.suptitle(f'{title_prefix}\nFeature Importance by Model (195 OTUs)',
                 fontsize=16, weight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(filename, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {filename}")


def draw_top50_grouped_circular(df_in, filename):
    """Top 50 OTU的综合环形图，4个模型用分组条形显示"""
    df_top = df_in.tail(50).reset_index(drop=True)
    n_bars = len(df_top)
    n_models = len(models)

    fig = plt.figure(figsize=(16, 16))
    ax = fig.add_subplot(111, polar=True)

    # 每组4个条形
    total_slots = n_bars * n_models
    theta_all = np.linspace(0, 2 * np.pi, total_slots, endpoint=False)
    width = 2 * np.pi / total_slots * 0.85

    for m_idx, model in enumerate(models):
        indices = np.arange(m_idx, total_slots, n_models)
        theta_model = theta_all[indices]
        radii = df_top[model].values

        ax.bar(theta_model, radii, width=width, bottom=0.02,
               color=colors_models[model], edgecolor='white', linewidth=0.3,
               alpha=0.85, label=model)

    # 给每个OTU组加一个标签（放在每组中间位置，即index 1或2）
    label_indices = np.arange(1, total_slots, n_models)
    for i in range(n_bars):
        angle = theta_all[label_indices[i]]
        max_r = df_top[models].iloc[i].max() + 0.1
        label = df_top['Label'].iloc[i]
        if len(label) > 20:
            label = label[:17] + '...'

        angle_deg = np.degrees(angle)
        if 90 < angle_deg <= 270:
            ha = 'right'
            angle_deg += 180
        else:
            ha = 'left'

        ax.text(angle, max_r, label,
                ha=ha, va='center', fontsize=7,
                rotation=angle_deg - 90,
                rotation_mode='anchor',
                color='#111111', weight='bold')

    ax.set_xticks([])
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=9, color='gray')
    ax.set_ylim(0, 1.45)
    ax.spines['polar'].set_visible(False)

    # 图例
    legend_elements = [Patch(facecolor=colors_models[m], edgecolor='white', label=m) for m in models]
    ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.25, 1.1), fontsize=11)

    plt.title('Top 50 OTUs Feature Importance\n(Grouped by 4 Models)', fontsize=16, pad=30, weight='bold')
    plt.savefig(filename, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {filename}")


if __name__ == '__main__':
    n_otus = len(df_sorted)
    print(f"Total OTUs after cleaning: {n_otus}")

    # 1. 综合重要性 - 全部OTU（平均重要性）
    draw_circular_bar_single(
        df_sorted, 'MeanImportance',
        'OTU综合重要性环形图 (195 OTUs, 4模型平均值)',
        'circular_bar_all_mean.svg',
        max_bars=n_otus, label_top_n=50
    )

    # 2. 综合重要性 - Top 100
    draw_circular_bar_single(
        df_sorted, 'MeanImportance',
        'OTU综合重要性环形图 (Top 100 OTUs, 4模型平均值)',
        'circular_bar_top100_mean.svg',
        max_bars=min(100, n_otus), label_top_n=50
    )

    # 3. 2x2分面 - 每个模型一个环（全部OTU）
    draw_circular_bar_facet(
        df_sorted,
        'OTU特征重要性',
        'circular_bar_facet_2x2_all.svg'
    )

    # 4. Top 50 分组环形图（4模型并排）
    draw_top50_grouped_circular(
        df_sorted,
        'circular_bar_top50_grouped.svg'
    )

    print("\n全部完成！")
