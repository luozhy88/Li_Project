import numpy as np

n_bars = 50
n_models = 4
total_slots = n_bars * n_models
theta_all = np.linspace(0, 2 * np.pi, total_slots, endpoint=False)
delta = 2 * np.pi / total_slots
width = delta * 0.82

# 检查第0组第一根柱子（红色，DRF）和前方分隔线的距离
sep_angle = theta_all[0] - (delta - width) / 2
sep_angle_norm = sep_angle % (2 * np.pi)
red_start = theta_all[0]
red_end = theta_all[0] + width
prev_purple_end = theta_all[199] + width

print(f'当前设置: width_ratio=0.82')
print(f'  每组间隙角度 = {(delta - width) * 180 / np.pi:.3f}°')
print(f'  分隔线到红色柱子边缘距离 = {((red_start - sep_angle_norm) % (2*np.pi)) * 180 / np.pi:.3f}°')
print()

# 如果改为 width_ratio=0.65
width2 = delta * 0.65
print(f'如果改为 width_ratio=0.65:')
print(f'  每组间隙角度 = {(delta - width2) * 180 / np.pi:.3f}°')
sep2 = theta_all[0] - (delta - width2) / 2
sep2_norm = sep2 % (2 * np.pi)
print(f'  分隔线到红色柱子边缘距离 = {((red_start - sep2_norm) % (2*np.pi)) * 180 / np.pi:.3f}°')
