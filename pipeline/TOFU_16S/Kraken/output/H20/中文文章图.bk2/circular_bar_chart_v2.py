import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, FancyArrowPatch
import os

# 字体
plt.rcParams['font.family'] = ['Noto Sans CJK JP', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 读取数据
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


def add_labels_with_leaders(ax, theta, radii, labels, label_indices, 
                            base_offset=0.12, min_arc_gap_deg=4.0, fontsize=7):
    """
    用引线方式添加标签，带简单的碰撞检测：
    - 标签按角度排序
    - 相邻标签如果角度太近，就交替放到不同的半径层
    - 用细灰线连接条形顶部和标签
    """
    if len(label_indices) == 0:
        return

    # 提取要标签的数据，按角度排序
    items = []
    for idx in label_indices:
        angle = theta[idx]
        radius = radii[idx]
        label = labels[idx]
        if len(label) > 22:
            label = label[:19] + '...'
        angle_deg = (np.degrees(angle) % 360)
        items.append({
            'idx': idx,
            'angle': angle,
            'angle_deg': angle_deg,
            'radius': radius,
            'label': label
        })
    
    # 按角度排序（顺时针）
    items.sort(key=lambda x: x['angle_deg'])
    
    # 计算每个标签应放的层级（简单的碰撞避免）
    # layer 0 = base_offset, layer 1 = base_offset + step, ...
    step = 0.10
    max_layer = 0
    for i, it in enumerate(items):
        if i == 0:
            it['layer'] = 0
            continue
        prev = items[i-1]
        arc_gap = abs(it['angle_deg'] - prev['angle_deg'])
        if arc_gap > 180:
            arc_gap = 360 - arc_gap
        if arc_gap < min_arc_gap_deg:
            # 跟上一个太近，放到下一层
            it['layer'] = prev['layer'] + 1
        else:
            it['layer'] = 0
        max_layer = max(max_layer, it['layer'])
    
    # 绘制标签和引线
    for it in items:
        angle = it['angle']
        radius = it['radius']
        label = it['label']
        layer = it['layer']
        label_r = 1.05 + base_offset + layer * step
        
        # 引线：从条形顶部到标签位置
        ax.plot([angle, angle], [radius + 0.02, label_r - 0.02], 
                color='gray', linewidth=0.5, alpha=0.6, zorder=2)
        
        angle_deg = np.degrees(angle)
        # 对齐：右半圆(90~270)文字朝内，左半圆朝外
        if 90 < angle_deg <= 270:
            ha = 'right'
            text_angle = angle_deg + 180
        else:
            ha = 'left'
            text_angle = angle_deg
        
        ax.text(angle, label_r, label,
                ha=ha, va='center', fontsize=fontsize,
                rotation=text_angle - 90,
                rotation_mode='anchor',
                color='#222222', weight='medium',
                zorder=3)


def draw_circular_bar_single(df_in, value_col, title, filename, 
                             max_bars=194, label_top_n=30):
    df_plot = df_in.tail(max_bars).reset_index(drop=True)
    n_bars = len(df_plot)
    
    # 增大画布，给外层标签留空间
    fig = plt.figure(figsize=(16, 16))
    ax = fig.add_subplot(111, polar=True)
    
    theta = np.linspace(0, 2 * np.pi, n_bars, endpoint=False)
    width = 2 * np.pi / n_bars * 0.82
    radii = df_plot[value_col].values
    
    norm = plt.Normalize(radii.min(), radii.max())
    cmap = plt.cm.viridis
    bar_colors = cmap(norm(radii))
    
    ax.bar(theta, radii, width=width, bottom=0.03, color=bar_colors, 
           edgecolor='white', linewidth=0.3, zorder=1)
    
    # 标签索引：只取最重要的
    label_indices = list(range(max(0, n_bars - label_top_n), n_bars))
    
    add_labels_with_leaders(ax, theta, radii, df_plot['Label'].values, 
                            label_indices, base_offset=0.08, fontsize=8)
    
    ax.set_xticks([])
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=9, color='gray')
    # 根据最大层数调整ylim
    ax.set_ylim(0, 1.55)
    ax.spines['polar'].set_visible(False)
    ax.grid(True, alpha=0.3)
    
    # 颜色条
    cax = fig.add_axes([0.88, 0.12, 0.025, 0.18])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label('Importance', fontsize=10)
    
    plt.title(title, fontsize=18, pad=35, weight='bold')
    plt.savefig(filename, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {filename}")


def draw_circular_bar_facet(df_in, title_prefix, filename):
    """2x2分面，增大间距，减少标签"""
    fig = plt.figure(figsize=(20, 20))
    
    n_bars = len(df_in)
    theta = np.linspace(0, 2 * np.pi, n_bars, endpoint=False)
    width = 2 * np.pi / n_bars * 0.88
    
    for idx, model in enumerate(models):
        # 使用gridspec手动布局，增加间距
        ax = fig.add_subplot(2, 2, idx+1, polar=True)
        radii = df_in[model].values
        
        norm = plt.Normalize(0, 1.0)
        cmap = plt.cm.plasma
        bar_colors = cmap(norm(radii))
        
        ax.bar(theta, radii, width=width, bottom=0.02, color=bar_colors, 
               edgecolor='white', linewidth=0.15, zorder=1)
        
        # 每个模型只标Top 20，用引线
        top_n = 20
        top_indices = np.argsort(radii)[-top_n:].tolist()
        
        add_labels_with_leaders(ax, theta, radii, df_in['Label'].values,
                                top_indices, base_offset=0.06, min_arc_gap_deg=5.0, fontsize=7)
        
        ax.set_xticks([])
        ax.set_yticks([0.25, 0.5, 0.75, 1.0])
        ax.set_yticklabels(['0.25', '0.5', '0.75', '1.0'], fontsize=8, color='gray')
        ax.set_ylim(0, 1.45)
        ax.spines['polar'].set_visible(False)
        ax.grid(True, alpha=0.25)
        ax.set_title(model, fontsize=15, pad=18, weight='bold', color=colors_models[model])
    
    plt.suptitle(f'{title_prefix}\nFeature Importance by Model (194 OTUs)',
                 fontsize=18, weight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96], h_pad=0.15, w_pad=0.15)
    plt.savefig(filename, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {filename}")


def draw_top50_grouped_circular(df_in, filename):
    """Top 50分组环形图，优化标签避免重叠"""
    df_top = df_in.tail(50).reset_index(drop=True)
    n_bars = len(df_top)
    n_models = len(models)
    
    fig = plt.figure(figsize=(18, 18))
    ax = fig.add_subplot(111, polar=True)
    
    total_slots = n_bars * n_models
    theta_all = np.linspace(0, 2 * np.pi, total_slots, endpoint=False)
    width = 2 * np.pi / total_slots * 0.82
    
    for m_idx, model in enumerate(models):
        indices = np.arange(m_idx, total_slots, n_models)
        theta_model = theta_all[indices]
        radii = df_top[model].values
        ax.bar(theta_model, radii, width=width, bottom=0.02,
               color=colors_models[model], edgecolor='white', linewidth=0.25,
               alpha=0.85, label=model, zorder=1)
    
    # 50个标签，每组一个，放在组中间（索引2）
    label_indices = np.arange(2, total_slots, n_models).tolist()
    # 计算每组的max radius
    max_radii = [df_top[models].iloc[i].max() for i in range(n_bars)]
    # 构建和 theta_all 等长的 radii 和 labels 数组
    radii_full = np.zeros(total_slots)
    labels_full = np.empty(total_slots, dtype=object)
    for i in range(n_bars):
        for j in range(n_models):
            slot = i * n_models + j
            radii_full[slot] = max_radii[i]
            labels_full[slot] = df_top['Label'].iloc[i]
    
    add_labels_with_leaders(ax, theta_all, radii_full, labels_full,
                            label_indices, base_offset=0.08, min_arc_gap_deg=6.0, fontsize=8)
    
    ax.set_xticks([])
    ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
    ax.set_yticklabels(['0.2', '0.4', '0.6', '0.8', '1.0'], fontsize=9, color='gray')
    ax.set_ylim(0, 1.55)
    ax.spines['polar'].set_visible(False)
    ax.grid(True, alpha=0.3)
    
    legend_elements = [Patch(facecolor=colors_models[m], edgecolor='white', label=m) for m in models]
    ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.22, 1.08), fontsize=12)
    
    plt.title('Top 50 OTUs Feature Importance\n(Grouped by 4 Models)', 
              fontsize=18, pad=35, weight='bold')
    plt.savefig(filename, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {filename}")


if __name__ == '__main__':
    n_otus = len(df_sorted)
    print(f"Total OTUs after cleaning: {n_otus}")
    
    draw_circular_bar_single(
        df_sorted, 'MeanImportance',
        'OTU综合重要性环形图 (194 OTUs, 4模型平均值)',
        'circular_bar_all_mean_v2.svg',
        max_bars=n_otus, label_top_n=30
    )
    
    draw_circular_bar_single(
        df_sorted, 'MeanImportance',
        'OTU综合重要性环形图 (Top 100 OTUs, 4模型平均值)',
        'circular_bar_top100_mean_v2.svg',
        max_bars=min(100, n_otus), label_top_n=30
    )
    
    draw_circular_bar_facet(
        df_sorted,
        'OTU特征重要性',
        'circular_bar_facet_2x2_all_v2.svg'
    )
    
    draw_top50_grouped_circular(
        df_sorted,
        'circular_bar_top50_grouped_v2.svg'
    )
    
    print("\n全部完成！")
