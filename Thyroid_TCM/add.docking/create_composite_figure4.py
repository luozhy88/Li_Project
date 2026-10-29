#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Composite docking figure for the main text: 2 rows x 2 columns.
Row A: 6UNI (human CYP3A4)   Row B: 9AYG (human CaV3.2/CACNA1H)
Col 1: local interaction     Col 2: global receptor view
The bacterial 4P5A entry is excluded per the manuscript.
"""

from PIL import Image, ImageDraw, ImageFont
import os

local_images = [
    'input/docking_6UNI_CID222284_20260423_132950/output_local.png',
    'input/docking_9AYG_CID222284_20260423_132658/output_local.png',
]
global_images = [
    'input/docking_6UNI_CID222284_20260423_132950/output_global.png',
    'input/docking_9AYG_CID222284_20260423_132658/output_global.png',
]
labels = ['A', 'B']
protein_names = ['6UNI', '9AYG']

thumb_w, thumb_h = 850, 638
label_w = 70
spacing_x = 25
spacing_y = 25
padding = 30

canvas_w = padding + label_w + thumb_w + spacing_x + thumb_w + padding
canvas_h = padding + (thumb_h + spacing_y) * 2 - spacing_y + padding

canvas = Image.new('RGB', (canvas_w, canvas_h), (255, 255, 255))
draw = ImageDraw.Draw(canvas)

label_font = None
title_font = None
for fp in ['/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf',
           '/System/Library/Fonts/Helvetica.ttc']:
    if os.path.exists(fp):
        label_font = ImageFont.truetype(fp, 56)
        title_font = ImageFont.truetype(fp, 32)
        break
if label_font is None:
    label_font = ImageFont.load_default()
    title_font = ImageFont.load_default()

for row in range(2):
    y_top = padding + row * (thumb_h + spacing_y)

    label = labels[row]
    bbox = draw.textbbox((0, 0), label, font=label_font)
    text_w = bbox[2] - bbox[0]
    text_h = bbox[3] - bbox[1]
    label_x = padding + (label_w - text_w) // 2
    label_y = y_top + (thumb_h - text_h) // 2
    draw.text((label_x, label_y), label, fill=(0, 0, 0), font=label_font)

    pname = protein_names[row]
    bbox2 = draw.textbbox((0, 0), pname, font=title_font)
    pw = bbox2[2] - bbox2[0]
    draw.text((padding + (label_w - pw) // 2, label_y + text_h + 10), pname, fill=(80, 80, 80), font=title_font)

    x1 = padding + label_w
    img_local = Image.open(local_images[row]).resize((thumb_w, thumb_h), Image.LANCZOS)
    canvas.paste(img_local, (x1, y_top))

    x2 = x1 + thumb_w + spacing_x
    img_global = Image.open(global_images[row]).resize((thumb_w, thumb_h), Image.LANCZOS)
    canvas.paste(img_global, (x2, y_top))

output_path = '../edit.Paper/output/20261024_03_support/20261024_Figure4.png'
canvas.save(output_path, dpi=(300, 300))
print(f'Composite image saved to: {output_path}')
print(f'Image size: {canvas.size}')
