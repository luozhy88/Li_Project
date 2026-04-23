#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create a composite image: 3 rows × 2 columns
Row 1: 4P5A (A)  Row 2: 6UNI (B)  Row 3: 9AYG (C)
Col 1: local interaction   Col 2: global docking
"""

from PIL import Image, ImageDraw, ImageFont
import os

# Config
local_images = [
    'input/docking_4P5A_CID222284_20260423_133247/output_local.png',
    'input/docking_6UNI_CID222284_20260423_132950/output_local.png',
    'input/docking_9AYG_CID222284_20260423_132658/output_local.png',
]
global_images = [
    'input/docking_4P5A_CID222284_20260423_133247/output_global.png',
    'input/docking_6UNI_CID222284_20260423_132950/output_global.png',
    'input/docking_9AYG_CID222284_20260423_132658/output_global.png',
]
labels = ['A', 'B', 'C']
protein_names = ['4P5A', '6UNI', '9AYG']

# Target size for each thumbnail (width, height)
thumb_w, thumb_h = 850, 638  # roughly 3:4 aspect from 2400x1800

label_w = 70      # width for row label (A, B, C)
spacing_x = 25    # horizontal gap between columns
spacing_y = 25    # vertical gap between rows
padding = 30      # outer padding

# Calculate canvas size
canvas_w = padding + label_w + thumb_w + spacing_x + thumb_w + padding
canvas_h = padding + (thumb_h + spacing_y) * 3 - spacing_y + padding

canvas = Image.new('RGB', (canvas_w, canvas_h), (255, 255, 255))
draw = ImageDraw.Draw(canvas)

# Try to load a font
try:
    # Try system fonts
    font_paths = [
        '/System/Library/Fonts/Helvetica.ttc',
        '/System/Library/Fonts/Arial.ttf',
        '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf',
    ]
    label_font = None
    title_font = None
    for fp in font_paths:
        if os.path.exists(fp):
            label_font = ImageFont.truetype(fp, 56)
            title_font = ImageFont.truetype(fp, 32)
            break
    if label_font is None:
        label_font = ImageFont.load_default()
        title_font = ImageFont.load_default()
except Exception:
    label_font = ImageFont.load_default()
    title_font = ImageFont.load_default()

for row in range(3):
    y_top = padding + row * (thumb_h + spacing_y)
    
    # Draw row label (A, B, C) on the left
    label = labels[row]
    bbox = draw.textbbox((0, 0), label, font=label_font)
    text_w = bbox[2] - bbox[0]
    text_h = bbox[3] - bbox[1]
    label_x = padding + (label_w - text_w) // 2
    label_y = y_top + (thumb_h - text_h) // 2
    draw.text((label_x, label_y), label, fill=(0, 0, 0), font=label_font)
    
    # Draw protein name under the label, smaller
    pname = protein_names[row]
    bbox2 = draw.textbbox((0, 0), pname, font=title_font)
    pw = bbox2[2] - bbox2[0]
    ph = bbox2[3] - bbox2[1]
    draw.text((padding + (label_w - pw) // 2, label_y + text_h + 10), pname, fill=(80, 80, 80), font=title_font)
    
    # Column 1: local interaction
    x1 = padding + label_w
    img_local = Image.open(local_images[row])
    img_local = img_local.resize((thumb_w, thumb_h), Image.LANCZOS)
    canvas.paste(img_local, (x1, y_top))
    
    # Column 2: global docking
    x2 = x1 + thumb_w + spacing_x
    img_global = Image.open(global_images[row])
    img_global = img_global.resize((thumb_w, thumb_h), Image.LANCZOS)
    canvas.paste(img_global, (x2, y_top))

# Save
output_path = 'output/figure_composite_docking.png'
canvas.save(output_path, dpi=(300, 300))
print(f'Composite image saved to: {output_path}')
print(f'Image size: {canvas.size}')
