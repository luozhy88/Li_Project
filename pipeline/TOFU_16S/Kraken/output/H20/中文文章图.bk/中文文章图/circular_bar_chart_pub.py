import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

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

model_labels = {
    'DRF': '分布式随机森林',
    'GBM': '梯度提升机',
    'GLM': '广义线性模型',
    'DeepLearning': '深度学习'
}

DPI = 600


def save_both(fig, filename_base):
    """同时保存 SVG 和 PDF"""
    for ext in ['.svg', '.pdf']:
        fname = filename_base + ext
        fig.savefig(fname, dpi=DPI, bbox_inches='tight',
                    facecolor='white', pad_inches=0.08)
        print(f"Saved: {fname}")


def draw_circular_bar_single(df_in, value_col, title, filename_base,
                             max_bars=194, label_top_n=20):
    df_plot = df_in.tail(max_bars).reset_index(drop=True)
    n_bars = len(df_plot)

    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(111, polar=True)

    theta = np.linspace(0, 2 * np.pi, n_bars, endpoint=False)
    width = 2 * np.pi / n_bars * 0.85
    radii = df_plot[value_col].values

    norm = plt.Normalize(radii.min(), radii.max())
    cmap = plt.cm.viridis
    bar_colors = cmap(norm(radii))

    ax.bar(theta, radii, width=width, bottom=0.02, color=bar_colors,
           edgecolor='white', linewidth=0.25)

    label_indices = list(range(max(0, n_bars - label_top_n), n_bars))

    for idx in label_indices:
        angle = theta[idx]
        radius = radii[idx] + 0.04
        label = df_plot['Label'].iloc[idx]
        if len(label) > 20:
            label = label[:17] + '...'

        angle_deg = np.degrees(angle)
        if 90 < angle_deg <= 270:
            ha = 'right'
            text_angle = angle_deg + 180
        else:
            ha = 'left'
            text_angle = angle_deg

        ax.text(angle, radius, label,
                ha=ha, va='center', fontsize=8,
                rotation=text_angle - 90,
                rotation_mode='anchor',
                color='#333333')

    ax.set_xticks([])
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=10, color='gray')
    ax.set_ylim(0, 1.12)
    ax.spines['polar'].set_visible(False)

    # 颜色条紧贴图
    cax = fig.add_axes([0.80, 0.28, 0.025, 0.12])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, cax=cax, label='Importance')

    plt.title(title, fontsize=16, pad=12, weight='bold')
    save_both(fig, filename_base)
    plt.close()


def draw_circular_bar_facet(df_in, title_prefix, filename_base):
    fig = plt.figure(figsize=(16, 16))

    n_bars = len(df_in)
    theta = np.linspace(0, 2 * np.pi, n_bars, endpoint=False)
    width = 2 * np.pi / n_bars * 0.9

    for idx, model in enumerate(models):
        ax = fig.add_subplot(2, 2, idx + 1, polar=True)
        radii = df_in[model].values

        norm = plt.Normalize(0, 1.0)
        cmap = plt.cm.plasma
        bar_colors = cmap(norm(radii))

        ax.bar(theta, radii, width=width, bottom=0.015,
               color=bar_colors, edgecolor='white', linewidth=0.15)

        top_n = 15
        top_indices = np.argsort(radii)[-top_n:]

        for ti in top_indices:
            angle = theta[ti]
            radius = radii[ti] + 0.04
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
                    ha=ha, va='center', fontsize=6,
                    rotation=text_angle - 90,
                    rotation_mode='anchor',
                    color='#222222')

        ax.set_xticks([])
        ax.set_yticks([0.25, 0.5, 0.75, 1.0])
        ax.set_yticklabels(['0.25', '0.5', '0.75', '1.0'], fontsize=8, color='gray')
        ax.set_ylim(0, 1.12)
        ax.spines['polar'].set_visible(False)
        ax.set_title(model_labels[model], fontsize=13, pad=8, weight='bold', color=colors_models[model])

    plt.suptitle(f'{title_prefix}\n模型特征重要性 (194 OTUs)',
                 fontsize=16, weight='bold', y=0.945)
    plt.tight_layout(rect=[0, 0, 1, 0.94], h_pad=0.08, w_pad=0.08)
    save_both(fig, filename_base)
    plt.close()


def draw_top50_grouped_circular(df_in, filename_base):
    df_top = df_in.tail(50).reset_index(drop=True)
    n_bars = len(df_top)
    n_models = len(models)

    fig = plt.figure(figsize=(13, 13))
    ax = fig.add_subplot(111, polar=True)

    total_slots = n_bars * n_models
    theta_all = np.linspace(0, 2 * np.pi, total_slots, endpoint=False)
    width = 2 * np.pi / total_slots * 0.85

    for m_idx, model in enumerate(models):
        indices = np.arange(m_idx, total_slots, n_models)
        theta_model = theta_all[indices]
        radii = df_top[model].values
        ax.bar(theta_model, radii, width=width, bottom=0.015,
               color=colors_models[model], edgecolor='white', linewidth=0.25,
               alpha=0.85, label=model)

    show_every = 2
    label_indices = np.arange(1, total_slots, n_models * show_every)

    for li in label_indices:
        group_idx = li // n_models
        angle = theta_all[li]
        max_r = df_top[models].iloc[group_idx].max() + 0.06
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
                ha=ha, va='center', fontsize=8,
                rotation=text_angle - 90,
                rotation_mode='anchor',
                color='#111111', weight='bold')

    ax.set_xticks([])
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=10, color='gray')
    ax.set_ylim(0, 1.22)
    ax.spines['polar'].set_visible(False)

    legend_elements = [Patch(facecolor=colors_models[m], edgecolor='white', label=model_labels[m]) for m in models]
    ax.legend(handles=legend_elements, loc='upper right',
              bbox_to_anchor=(1.08, 1.02), fontsize=11, frameon=False)

    plt.title('',
              fontsize=16, pad=12, weight='bold')
    save_both(fig, filename_base)
    plt.close()


if __name__ == '__main__':
    n_otus = len(df_sorted)
    print(f"Total OTUs after cleaning: {n_otus}")

    draw_circular_bar_single(
        df_sorted, 'MeanImportance',
        'OTU综合重要性环形图 (194 OTUs, 4模型平均值)',
        'circular_bar_all_mean_pub',
        max_bars=n_otus, label_top_n=20
    )

    draw_circular_bar_single(
        df_sorted, 'MeanImportance',
        'OTU综合重要性环形图 (Top 100 OTUs, 4模型平均值)',
        'circular_bar_top100_mean_pub',
        max_bars=min(100, n_otus), label_top_n=20
    )

    draw_circular_bar_facet(
        df_sorted,
        'OTU特征重要性',
        'circular_bar_facet_2x2_all_pub'
    )

    draw_top50_grouped_circular(
        df_sorted,
        'circular_bar_top50_grouped_pub'
    )

    print("\n全部完成！")
