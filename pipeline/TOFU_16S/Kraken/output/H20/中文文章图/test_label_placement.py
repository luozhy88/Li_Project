import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

fig, axes = plt.subplots(2, 2, figsize=(12, 12), subplot_kw=dict(polar=True))

angles = np.radians([0, 45, 90, 135, 180, 225, 270, 315])
labels = ['0°', '45°', '90°', '135°', '180°', '225°', '270°', '315°']

for idx, (title, inner) in enumerate([('外侧径向', False), ('内侧径向', True), 
                                       ('外侧居中', False), ('内侧居中', True)]):
    ax = axes[idx // 2, idx % 2]
    ax.set_ylim(0, 1.2)
    
    for angle_deg, label in zip(np.degrees(angles), labels):
        angle_rad = np.radians(angle_deg)
        
        if inner:
            radius = 0.5
        else:
            radius = 1.0
        
        if '径向' in title:
            # 径向排列
            if 90 < angle_deg <= 270:
                rot = angle_deg + 180 if not inner else angle_deg
                ha = 'right' if not inner else 'left'
            else:
                rot = angle_deg if not inner else angle_deg + 180
                ha = 'left' if not inner else 'right'
        else:
            # 居中径向
            rot = angle_deg if angle_deg <= 90 or angle_deg > 270 else angle_deg + 180
            ha = 'center'
        
        ax.text(angle_rad, radius, label, ha=ha, va='center', 
                rotation=rot, fontsize=10, color='black')
        ax.plot([angle_rad, angle_rad], [0, 1.1], 'k-', linewidth=0.5, alpha=0.3)
    
    ax.set_title(title, fontsize=12)

plt.tight_layout()
plt.savefig('test_label_placement.pdf', dpi=150)
plt.close()
print('Saved test_label_placement.pdf')
