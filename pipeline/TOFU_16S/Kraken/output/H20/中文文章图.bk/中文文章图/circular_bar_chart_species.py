import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

plt.rcParams['font.family'] = ['Noto Sans CJK JP', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 读取种级别数据
df = pd.read_csv('../SRP151288_HC_vs_TC.Species/H20_feature_importance.csv')
models = ['DRF', 'GBM', 'GLM']
df = df.dropna(subset=models, how='all').reset_index(drop=True)
df['MeanImportance'] = df[models].mean(axis=1)
df['MaxImportance'] = df[models].max(axis=1)
df_sorted = df.sort_values('MeanImportance', ascending=True).reset_index(drop=True)
# 提取种名：去掉 OTUxxx_s__ 前缀
df_sorted['Label'] = df_sorted['variable'].str.replace(r'^OTU\d+_s__', '', regex=True)

colors_models = {
    'DRF': '#E74C3C',
    'GBM': '#3498DB',
    'GLM': '#2ECC71'
}

DPI = 600


def save_both(fig, filename_base):
    for ext in ['.svg', '.pdf']:
        fname = filename_base + ext
        fig.savefig(fname, dpi=DPI, bbox_inches='tight',
                    facecolor='white', pad_inches=0.08)
        print(f"Saved: {fname}")


def draw_circular_bar_single(df_in, value_col, title, filename_base,
                             max_bars=None, label_top_n=None):
    if max_bars is None:
        max_bars = len(df_in)
    if label_top_n is None:
        label_top_n = len(df_in)

    df_plot = df_in.tail(max_bars).reset_index(drop=True)
    n_bars = len(df_plot)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, polar=True)

    theta = np.linspace(0, 2 * np.pi, n_bars, endpoint=False)
    width = 2 * np.pi / n_bars * 0.82
    radii = df_plot[value_col].values

    norm = plt.Normalize(radii.min(), radii.max())
    cmap = plt.cm.viridis
    bar_colors = cmap(norm(radii))

    ax.bar(theta, radii, width=width, bottom=0.02, color=bar_colors,
           edgecolor='white', linewidth=0.4)

    label_indices = list(range(max(0, n_bars - label_top_n), n_bars))

    for idx in label_indices:
        angle = theta[idx]
        radius = radii[idx] + 0.05
        label = df_plot['Label'].iloc[idx]
        if len(label) > 25:
            label = label[:22] + '...'

        angle_deg = np.degrees(angle)
        if 90 < angle_deg <= 270:
            ha = 'right'
            text_angle = angle_deg + 180
        else:
            ha = 'left'
            text_angle = angle_deg

        ax.text(angle, radius, label,
                ha=ha, va='center', fontsize=9,
                rotation=text_angle - 90,
                rotation_mode='anchor',
                color='#333333')

    ax.set_xticks([])
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=10, color='gray')
    ax.set_ylim(0, 1.25)
    ax.spines['polar'].set_visible(False)

    cax = fig.add_axes([0.78, 0.28, 0.025, 0.14])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.colorbar(sm, cax=cax, label='Importance')

    plt.title(title, fontsize=15, pad=12, weight='bold')
    save_both(fig, filename_base)
    plt.close()


def draw_circular_bar_facet(df_in, title_prefix, filename_base):
    n_bars = len(df_in)
    theta = np.linspace(0, 2 * np.pi, n_bars, endpoint=False)
    width = 2 * np.pi / n_bars * 0.88

    fig = plt.figure(figsize=(16, 6))

    for idx, model in enumerate(models):
        ax = fig.add_subplot(1, 3, idx + 1, polar=True)
        radii = df_in[model].values

        norm = plt.Normalize(0, 1.0)
        cmap = plt.cm.plasma
        bar_colors = cmap(norm(radii))

        ax.bar(theta, radii, width=width, bottom=0.02,
               color=bar_colors, edgecolor='white', linewidth=0.3)

        # 22个OTU可以全部标标签
        for ti in range(n_bars):
            angle = theta[ti]
            radius = radii[ti] + 0.05
            label = df_in['Label'].iloc[ti]
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
                    color='#222222')

        ax.set_xticks([])
        ax.set_yticks([0.25, 0.5, 0.75, 1.0])
        ax.set_yticklabels(['0.25', '0.5', '0.75', '1.0'], fontsize=8, color='gray')
        ax.set_ylim(0, 1.25)
        ax.spines['polar'].set_visible(False)
        ax.set_title(model, fontsize=13, pad=10, weight='bold', color=colors_models[model])

    plt.suptitle(f'{title_prefix}\nFeature Importance by Model ({n_bars} OTUs)',
                 fontsize=15, weight='bold', y=1.02)
    plt.tight_layout(rect=[0, 0, 1, 0.98], h_pad=0.1, w_pad=0.1)
    save_both(fig, filename_base)
    plt.close()


def draw_grouped_circular(df_in, filename_base):
    n_bars = len(df_in)
    n_models = len(models)

    fig = plt.figure(figsize=(11, 11))
    ax = fig.add_subplot(111, polar=True)

    total_slots = n_bars * n_models
    theta_all = np.linspace(0, 2 * np.pi, total_slots, endpoint=False)
    width = 2 * np.pi / total_slots * 0.82

    for m_idx, model in enumerate(models):
        indices = np.arange(m_idx, total_slots, n_models)
        theta_model = theta_all[indices]
        radii = df_in[model].values
        ax.bar(theta_model, radii, width=width, bottom=0.015,
               color=colors_models[model], edgecolor='white', linewidth=0.3,
               alpha=0.85, label=model)

    # 每组一个标签，放在中间位置（index 1）
    label_indices = np.arange(1, total_slots, n_models)
    radii_full = np.zeros(total_slots)
    labels_full = np.empty(total_slots, dtype=object)
    for i in range(n_bars):
        for j in range(n_models):
            slot = i * n_models + j
            radii_full[slot] = df_in[models].iloc[i].max()
            labels_full[slot] = df_in['Label'].iloc[i]

    for li in label_indices:
        angle = theta_all[li]
        radius = radii_full[li] + 0.07
        label = labels_full[li]
        if len(label) > 25:
            label = label[:22] + '...'

        angle_deg = np.degrees(angle)
        if 90 < angle_deg <= 270:
            ha = 'right'
            text_angle = angle_deg + 180
        else:
            ha = 'left'
            text_angle = angle_deg

        ax.text(angle, radius, label,
                ha=ha, va='center', fontsize=8.5,
                rotation=text_angle - 90,
                rotation_mode='anchor',
                color='#111111', weight='bold')

    ax.set_xticks([])
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=10, color='gray')
    ax.set_ylim(0, 1.28)
    ax.spines['polar'].set_visible(False)

    legend_elements = [Patch(facecolor=colors_models[m], edgecolor='white', label=m) for m in models]
    ax.legend(handles=legend_elements, loc='upper right',
              bbox_to_anchor=(1.08, 1.02), fontsize=11, frameon=False)

    plt.title(f'{n_bars} OTUs Feature Importance\n(Grouped by 3 Models)',
              fontsize=15, pad=12, weight='bold')
    save_both(fig, filename_base)
    plt.close()


if __name__ == '__main__':
    n_otus = len(df_sorted)
    print(f"Total species-level OTUs: {n_otus}")

    draw_circular_bar_single(
        df_sorted, 'MeanImportance',
        '种级别OTU综合重要性环形图 (3模型平均值)',
        'circular_bar_species_all_mean_pub',
        max_bars=n_otus, label_top_n=n_otus
    )

    draw_circular_bar_facet(
        df_sorted,
        '种级别OTU特征重要性',
        'circular_bar_species_facet_pub'
    )

    draw_grouped_circular(
        df_sorted,
        'circular_bar_species_grouped_pub'
    )

    print("\n全部完成！")
