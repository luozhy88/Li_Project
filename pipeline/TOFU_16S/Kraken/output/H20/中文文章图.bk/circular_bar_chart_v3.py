import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os

plt.rcParams['font.family'] = ['Noto Sans CJK JP', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

df = pd.read_csv('H20_feature_importance_data.csv')
models = ['DRF', 'GBM', 'GLM', 'DeepLearning']
df = df.dropna(subset=models, how='all').reset_index(drop=True)
df['MeanImportance'] = df[models].mean(axis=1)
df['MaxImportance'] = df[models].max(axis=1)
df_sorted = df.sort_values('MeanImportance', ascending=True).reset_index(drop=True)
df_sorted['Label'] = df_sorted['OTU'].str.replace(r'^OTU\d+_', '', regex=True)

colors_models = {
    'DRF': '#E74C3C',
    'GBM': '#3498DB',
    'GLM': '#2ECC71',
    'DeepLearning': '#9B59B6'
}


def draw_circular_bar_single(df_in, value_col, title, filename,
                             max_bars=194, label_top_n=20):
    df_plot = df_in.tail(max_bars).reset_index(drop=True)
    n_bars = len(df_plot)

    fig = plt.figure(figsize=(16, 16))
    ax = fig.add_subplot(111, polar=True)

    theta = np.linspace(0, 2 * np.pi, n_bars, endpoint=False)
    width = 2 * np.pi / n_bars * 0.85
    radii = df_plot[value_col].values

    norm = plt.Normalize(radii.min(), radii.max())
    cmap = plt.cm.viridis
    bar_colors = cmap(norm(radii))

    ax.bar(theta, radii, width=width, bottom=0.03, color=bar_colors,
           edgecolor='white', linewidth=0.3)

    # 只给尾部最重要的 label_top_n 个加标签
    label_indices = list(range(max(0, n_bars - label_top_n), n_bars))

    for idx in label_indices:
        angle = theta[idx]
        radius = radii[idx] + 0.06
        label = df_plot['Label'].iloc[idx]
        if len(label) > 20:
            label = label[:17] + '...'

        angle_deg = np.degrees(angle)
        # 右半圆标签放在条形右侧，文字朝外（顺着圆的切线）
        if 90 < angle_deg <= 270:
            ha = 'right'
            text_angle = angle_deg + 180
        else:
            ha = 'left'
            text_angle = angle_deg

        ax.text(angle, radius, label,
                ha=ha, va='center', fontsize=7.5,
                rotation=text_angle - 90,
                rotation_mode='anchor',
                color='#333333')

    ax.set_xticks([])
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=9, color='gray')
    ax.set_ylim(0, 1.18)
    ax.spines['polar'].set_visible(False)

    # 颜色条靠近图
    cax = fig.add_axes([0.82, 0.22, 0.025, 0.14])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, cax=cax, label='Importance')

    plt.title(title, fontsize=16, pad=18, weight='bold')
    plt.savefig(filename, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {filename}")


def draw_circular_bar_facet(df_in, title_prefix, filename):
    fig = plt.figure(figsize=(20, 20))

    n_bars = len(df_in)
    theta = np.linspace(0, 2 * np.pi, n_bars, endpoint=False)
    width = 2 * np.pi / n_bars * 0.9

    for idx, model in enumerate(models):
        ax = fig.add_subplot(2, 2, idx + 1, polar=True)
        radii = df_in[model].values

        norm = plt.Normalize(0, 1.0)
        cmap = plt.cm.plasma
        bar_colors = cmap(norm(radii))

        ax.bar(theta, radii, width=width, bottom=0.02,
               color=bar_colors, edgecolor='white', linewidth=0.2)

        # 每环只标 Top 15，减少重叠
        top_n = 15
        top_indices = np.argsort(radii)[-top_n:]

        for ti in top_indices:
            angle = theta[ti]
            radius = radii[ti] + 0.05
            label = df_in['Label'].iloc[ti]
            if len(label) > 18:
                label = label[:15] + '...'

            angle_deg = np.degrees(angle)
            if 90 < angle_deg <= 270:
                ha = 'right'
                text_angle = angle_deg + 180
            else:
                ha = 'left'
                text_angle = angle_deg

            ax.text(angle, radius, label,
                    ha=ha, va='center', fontsize=5.5,
                    rotation=text_angle - 90,
                    rotation_mode='anchor',
                    color='#222222')

        ax.set_xticks([])
        ax.set_yticks([0.25, 0.5, 0.75, 1.0])
        ax.set_yticklabels(['0.25', '0.5', '0.75', '1.0'], fontsize=7, color='gray')
        ax.set_ylim(0, 1.15)
        ax.spines['polar'].set_visible(False)
        ax.set_title(model, fontsize=14, pad=10, weight='bold', color=colors_models[model])

    plt.suptitle(f'{title_prefix}\nFeature Importance by Model (194 OTUs)',
                 fontsize=16, weight='bold', y=0.96)
    plt.tight_layout(rect=[0, 0, 1, 0.96], h_pad=0.12, w_pad=0.12)
    plt.savefig(filename, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {filename}")


def draw_top50_grouped_circular(df_in, filename):
    df_top = df_in.tail(50).reset_index(drop=True)
    n_bars = len(df_top)
    n_models = len(models)

    fig = plt.figure(figsize=(17, 17))
    ax = fig.add_subplot(111, polar=True)

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

    # 50个标签，每隔一个显示，进一步降低密度 -> 25个标签
    show_every = 2
    label_indices = np.arange(1, total_slots, n_models * show_every)

    for li in label_indices:
        group_idx = li // n_models
        angle = theta_all[li]
        max_r = df_top[models].iloc[group_idx].max() + 0.08
        label = df_top['Label'].iloc[group_idx]
        if len(label) > 20:
            label = label[:17] + '...'

        angle_deg = np.degrees(angle)
        if 90 < angle_deg <= 270:
            ha = 'right'
            text_angle = angle_deg + 180
        else:
            ha = 'left'
            text_angle = angle_deg

        ax.text(angle, max_r, label,
                ha=ha, va='center', fontsize=7.5,
                rotation=text_angle - 90,
                rotation_mode='anchor',
                color='#111111', weight='bold')

    ax.set_xticks([])
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=9, color='gray')
    ax.set_ylim(0, 1.28)
    ax.spines['polar'].set_visible(False)

    legend_elements = [Patch(facecolor=colors_models[m], edgecolor='white', label=m) for m in models]
    ax.legend(handles=legend_elements, loc='upper right',
              bbox_to_anchor=(1.12, 1.04), fontsize=11)

    plt.title('Top 50 OTUs Feature Importance\n(Grouped by 4 Models)',
              fontsize=16, pad=18, weight='bold')
    plt.savefig(filename, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {filename}")


if __name__ == '__main__':
    n_otus = len(df_sorted)
    print(f"Total OTUs after cleaning: {n_otus}")

    draw_circular_bar_single(
        df_sorted, 'MeanImportance',
        'OTU综合重要性环形图 (194 OTUs, 4模型平均值)',
        'circular_bar_all_mean_v3.svg',
        max_bars=n_otus, label_top_n=20
    )

    draw_circular_bar_single(
        df_sorted, 'MeanImportance',
        'OTU综合重要性环形图 (Top 100 OTUs, 4模型平均值)',
        'circular_bar_top100_mean_v3.svg',
        max_bars=min(100, n_otus), label_top_n=20
    )

    draw_circular_bar_facet(
        df_sorted,
        'OTU特征重要性',
        'circular_bar_facet_2x2_all_v3.svg'
    )

    draw_top50_grouped_circular(
        df_sorted,
        'circular_bar_top50_grouped_v3.svg'
    )

    print("\n全部完成！")
