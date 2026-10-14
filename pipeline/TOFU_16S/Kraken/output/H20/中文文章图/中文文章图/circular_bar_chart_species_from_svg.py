import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

plt.rcParams['font.family'] = ['Noto Sans CJK JP', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 读取从SVG解析出的种级数据
df = pd.read_csv('H20_feature_importance_species_from_svg.csv')
models = ['DRF', 'GBM', 'GLM', 'DeepLearning']
df = df.dropna(subset=models, how='all').reset_index(drop=True)
df['MeanImportance'] = df[models].mean(axis=1)
df['MaxImportance'] = df[models].max(axis=1)
df_sorted = df.sort_values('MeanImportance', ascending=True).reset_index(drop=True)
df_sorted['Label'] = df_sorted['OTU'].str.replace(r'^OTU\d+_s__', '', regex=True)

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
    for ext in ['.svg', '.pdf']:
        fname = filename_base + ext
        fig.savefig(fname, dpi=DPI, bbox_inches='tight',
                    facecolor='white', pad_inches=0.08)
        print(f"Saved: {fname}")


def draw_circular_bar_single(df_in, value_col, title, filename_base,
                             max_bars=None, label_top_n=25):
    if max_bars is None:
        max_bars = len(df_in)
    df_plot = df_in.tail(max_bars).reset_index(drop=True)
    n_bars = len(df_plot)

    fig = plt.figure(figsize=(18, 18))
    ax = fig.add_subplot(111, polar=True)

    theta = np.linspace(0, 2 * np.pi, n_bars, endpoint=False)
    width = 2 * np.pi / n_bars * 0.82
    radii = df_plot[value_col].values

    norm = plt.Normalize(radii.min(), radii.max())
    cmap = plt.cm.viridis
    bar_colors = cmap(norm(radii))

    ax.bar(theta, radii, width=width, bottom=0.0, color=bar_colors,
           edgecolor='white', linewidth=0.25)

    label_indices = list(range(max(0, n_bars - label_top_n), n_bars))

    for idx in label_indices:
        angle = theta[idx]
        radius = radii[idx] + 0.03
        label = df_plot['Label'].iloc[idx]
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
                color='#333333')

    ax.set_xticks([])
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=10, color='gray')
    ax.set_ylim(0, 1.05)
    ax.spines['polar'].set_visible(False)

    cax = fig.add_axes([0.85, 0.25, 0.022, 0.10])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, cax=cax, label='重要性')

    plt.title('', fontsize=16, pad=14, weight='bold')
    save_both(fig, filename_base)
    plt.close()


def draw_circular_bar_facet(df_in, title_prefix, filename_base):
    fig = plt.figure(figsize=(20, 20))

    n_bars = len(df_in)
    theta = np.linspace(0, 2 * np.pi, n_bars, endpoint=False)
    width = 2 * np.pi / n_bars * 0.88

    for idx, model in enumerate(models):
        ax = fig.add_subplot(2, 2, idx + 1, polar=True)
        radii = df_in[model].fillna(0).values

        norm = plt.Normalize(0, 1.0)
        cmap = plt.cm.plasma
        bar_colors = cmap(norm(radii))

        ax.bar(theta, radii, width=width, bottom=0.015,
               color=bar_colors, edgecolor='white', linewidth=0.15)

        top_n = 18
        top_indices = np.argsort(radii)[-top_n:]

        for ti in top_indices:
            angle = theta[ti]
            radius = radii[ti] + 0.04
            label = df_in['Label'].iloc[ti]
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
                    ha=ha, va='center', fontsize=5.5,
                    rotation=text_angle - 90,
                    rotation_mode='anchor',
                    color='#222222')

        ax.set_xticks([])
        ax.set_yticks([0.25, 0.5, 0.75, 1.0])
        ax.set_yticklabels(['0.25', '0.5', '0.75', '1.0'], fontsize=8, color='gray')
        ax.set_ylim(0, 1.18)
        ax.spines['polar'].set_visible(False)
        ax.set_title(model_labels[model], fontsize=13, pad=10, weight='bold', color=colors_models[model])

    plt.suptitle(f'{title_prefix}\n模型特征重要性 ({n_bars} OTUs)',
                 fontsize=16, weight='bold', y=0.96)
    plt.tight_layout(rect=[0, 0, 1, 0.95], h_pad=0.10, w_pad=0.10)
    save_both(fig, filename_base)
    plt.close()


def draw_grouped_circular(df_in, filename_base, top_n=50):
    df_top = df_in.tail(top_n).reset_index(drop=True)
    n_bars = len(df_top)
    n_models = len(models)

    fig = plt.figure(figsize=(16, 16))
    ax = fig.add_subplot(111, polar=True)

    total_slots = n_bars * n_models
    theta_all = np.linspace(0, 2 * np.pi, total_slots, endpoint=False)
    width = 2 * np.pi / total_slots * 0.82

    for m_idx, model in enumerate(models):
        indices = np.arange(m_idx, total_slots, n_models)
        theta_model = theta_all[indices]
        radii = df_top[model].fillna(0).values
        ax.bar(theta_model, radii, width=width, bottom=0.015,
               color=colors_models[model], edgecolor='white', linewidth=0.25,
               alpha=0.85, label=model)

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
        radius = radii_full[li] + 0.06
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
                ha=ha, va='center', fontsize=7.5,
                rotation=text_angle - 90,
                rotation_mode='anchor',
                color='#111111', weight='bold')

    ax.set_xticks([])
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=10, color='gray')
    ax.set_ylim(0, 1.28)
    ax.spines['polar'].set_visible(False)

    legend_elements = [Patch(facecolor=colors_models[m], edgecolor='white', label=model_labels[m]) for m in models]
    ax.legend(handles=legend_elements, loc='upper right',
              bbox_to_anchor=(1.04, 1.01), fontsize=10, frameon=False)

    plt.title('',
              fontsize=16, pad=14, weight='bold')
    save_both(fig, filename_base)
    plt.close()


if __name__ == '__main__':
    n_otus = len(df_sorted)
    print(f"Total species-level OTUs from SVG: {n_otus}")

    draw_circular_bar_single(
        df_sorted, 'MeanImportance',
        '种级别OTU综合重要性环形图 (351 OTUs, 4模型平均值)',
        'circular_bar_species_all_mean_svg',
        max_bars=n_otus, label_top_n=25
    )

    draw_circular_bar_single(
        df_sorted, 'MeanImportance',
        '种级别OTU综合重要性环形图 (Top 100, 4模型平均值)',
        'circular_bar_species_top100_mean_svg',
        max_bars=min(100, n_otus), label_top_n=25
    )

    draw_circular_bar_facet(
        df_sorted,
        '种级别OTU特征重要性',
        'circular_bar_species_facet_svg'
    )

    draw_grouped_circular(
        df_sorted,
        'circular_bar_species_top50_grouped_svg',
        top_n=50
    )

    print("\n全部完成！")
