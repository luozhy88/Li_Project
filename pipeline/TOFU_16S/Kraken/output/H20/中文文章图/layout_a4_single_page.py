import xml.etree.ElementTree as ET
import copy
import re

# Parse original SVG
svg_file = "H20_feature_importance_heatmap_plot属.svg"
tree = ET.parse(svg_file)
root = tree.getroot()

ET.register_namespace('', 'http://www.w3.org/2000/svg')
ET.register_namespace('xlink', 'http://www.w3.org/1999/xlink')

# Extract OTU labels
otu_texts = []
for elem in root.iter():
    if elem.tag == '{http://www.w3.org/2000/svg}text' and elem.get('text-anchor') == 'end' and elem.get('x') == '194.15':
        y = float(elem.get('y', 0))
        content = ''.join(elem.itertext())
        # Shorten label: remove OTUxxxx_ prefix
        short_content = re.sub(r'^OTU\d+_', '', content)
        otu_texts.append((y, elem, content, short_content))

otu_texts.sort(key=lambda x: x[0])
n_rows = len(otu_texts)
print(f"Total OTU rows: {n_rows}")

# Group into 3 columns
n_groups = 3
rows_per_group = n_rows // n_groups
remainder = n_rows % n_groups

groups = []
start = 0
for i in range(n_groups):
    size = rows_per_group + (1 if i < remainder else 0)
    groups.append(otu_texts[start:start+size])
    start += size

print(f"Groups: {[len(g) for g in groups]}")

# Collect heatmap rects
heatmap_rects = []
for elem in root.iter():
    if elem.tag == '{http://www.w3.org/2000/svg}rect' and elem.get('width') == '33.92':
        y = float(elem.get('y', 0))
        heatmap_rects.append((y, elem))

# Collect grid polylines
grid_polylines = []
for elem in root.iter():
    if elem.tag == '{http://www.w3.org/2000/svg}polyline':
        pts = elem.get('points', '')
        coords = re.findall(r'[\d.]+', pts)
        if len(coords) >= 2:
            y_val = float(coords[1])
            if not ('358' in pts or '375' in pts or '371' in pts or '361' in pts):
                grid_polylines.append((y_val, elem))

# Collect column labels
col_labels = []
for elem in root.iter():
    if elem.tag == '{http://www.w3.org/2000/svg}text' and 'transform' in elem.attrib and 'rotate' in elem.attrib.get('transform', ''):
        col_labels.append(elem)

# Collect colorbar elements
cbar_image = None
cbar_polylines = []
cbar_texts = []
title_text = None

for elem in root.iter():
    if elem.tag == '{http://www.w3.org/2000/svg}image':
        cbar_image = elem
    elif elem.tag == '{http://www.w3.org/2000/svg}text':
        content = ''.join(elem.itertext())
        x = elem.get('x')
        if x == '380.76' and content in ['0.00', '0.25', '0.50', '0.75', '1.00']:
            cbar_texts.append((content, elem))
        elif content == '特征重要性':
            title_text = elem
    elif elem.tag == '{http://www.w3.org/2000/svg}polyline':
        pts = elem.get('points', '')
        if '358' in pts or '375' in pts or '371' in pts or '361' in pts:
            cbar_polylines.append(elem)

order = {'0.00': 0, '0.25': 1, '0.50': 2, '0.75': 3, '1.00': 4}
cbar_texts.sort(key=lambda x: order.get(x[0], 99))

# A4 Landscape single page layout
A4_W, A4_H = 842, 595
MARGIN_L, MARGIN_R = 45, 45
MARGIN_T, MARGIN_B = 40, 35
usable_w = A4_W - MARGIN_L - MARGIN_R  # 752
usable_h = A4_H - MARGIN_T - MARGIN_B  # 520

n_cols = 3
col_gap = 15
col_w = (usable_w - col_gap * (n_cols - 1)) / n_cols  # ~237pt

# Original content width estimate after shortening labels
# Label area: original 194pt, after shortening ~120pt needed
# Heatmap area: 143pt
# Total: ~263pt
base_scale = col_w / 263.0
print(f"Column scale: {base_scale:.3f}")

row_height_orig = 9.1
row_height = row_height_orig * base_scale
max_rows = max(len(g) for g in groups)
content_h = max_rows * row_height
print(f"Row height: {row_height:.2f}pt, content height: {content_h:.1f}pt")

# Check if fits
bottom_margin_for_colnames = 25
total_h_needed = content_h + bottom_margin_for_colnames
print(f"Total height needed: {total_h_needed:.1f}pt, available: {usable_h}pt")

# Adjust if needed
if total_h_needed > usable_h:
    # Reduce scale to fit height
    scale_from_h = (usable_h - bottom_margin_for_colnames) / (max_rows * row_height_orig)
    base_scale = min(base_scale, scale_from_h)
    row_height = row_height_orig * base_scale
    content_h = max_rows * row_height
    total_h_needed = content_h + bottom_margin_for_colnames
    print(f"Adjusted scale: {base_scale:.3f}, total height: {total_h_needed:.1f}pt")

# Build single-page A4 landscape SVG
new_root = ET.Element('svg')
new_root.set('class', 'svglite')
new_root.set('width', f'{A4_W}pt')
new_root.set('height', f'{A4_H}pt')
new_root.set('viewBox', f'0 0 {A4_W} {A4_H}')

# Add defs with style
defs = ET.SubElement(new_root, 'defs')
style = ET.SubElement(defs, 'style')
style.set('type', 'text/css')
style.text = '''
    .svglite line, .svglite polyline, .svglite polygon, .svglite path, .svglite rect, .svglite circle {
      fill: none;
      stroke: #000000;
      stroke-linecap: round;
      stroke-linejoin: round;
      stroke-miterlimit: 10.00;
    }
    .svglite text {
      white-space: pre;
    }
  '''

# Background
bg = ET.SubElement(new_root, 'rect')
bg.set('width', '100%')
bg.set('height', '100%')
bg.set('style', 'stroke: none; fill: #FFFFFF;')

def shift_scale_shorten(elem, y_offset, x_offset, scale, new_text=None):
    """Apply y shift, x shift, uniform scale, and optionally replace text + remove textLength"""
    e = copy.deepcopy(elem)
    
    # Replace text content if provided
    if new_text is not None:
        # Clear existing text content
        for child in list(e):
            e.remove(child)
        if len(e) == 0 or (e.text is None and len(list(e)) == 0):
            e.text = new_text
        else:
            # More complex: replace text nodes
            e.text = new_text
            for child in list(e):
                e.remove(child)
        # Remove textLength since text changed
        if 'textLength' in e.attrib:
            del e.attrib['textLength']
    
    # Scale and shift y
    if e.get('y') is not None:
        try:
            old_y = float(e.get('y'))
            new_y = (old_y - y_offset) * scale
            e.set('y', f'{new_y:.4f}'.rstrip('0').rstrip('.'))
        except:
            pass
    
    # Scale and shift x
    if e.get('x') is not None:
        try:
            old_x = float(e.get('x'))
            new_x = old_x * scale + x_offset
            e.set('x', f'{new_x:.4f}'.rstrip('0').rstrip('.'))
        except:
            pass
    
    # Scale width and height
    for attr in ['width', 'height']:
        if e.get(attr) is not None:
            try:
                old_val = float(e.get(attr))
                new_val = old_val * scale
                e.set(attr, f'{new_val:.4f}'.rstrip('0').rstrip('.'))
            except:
                pass
    
    # Scale points in polyline
    if e.get('points') is not None:
        pts = e.get('points')
        coords = re.findall(r'([\d.]+),([\d.]+)', pts)
        new_pts = []
        for x_str, y_str in coords:
            new_x = float(x_str) * scale + x_offset
            new_y = (float(y_str) - y_offset) * scale
            new_pts.append(f"{new_x:.4f},{new_y:.4f}".rstrip('0').rstrip('.').replace(',.', ',0.').replace('.,', '.0,'))
        e.set('points', ' '.join(new_pts))
    
    # Scale transform
    if e.get('transform') is not None:
        trans = e.get('transform')
        def scale_translate(m):
            tx = float(m.group(1)) * scale + x_offset
            ty = (float(m.group(2)) - y_offset) * scale
            return f'translate({tx:.2f},{ty:.2f})'
        new_trans = re.sub(r'translate\(([\d.]+),\s*([\d.]+)\)', scale_translate, trans)
        e.set('transform', new_trans)
    
    # Scale textLength (if not removed)
    if e.get('textLength') is not None:
        try:
            tl = e.get('textLength').replace('px', '')
            new_tl = float(tl) * scale
            e.set('textLength', f'{new_tl:.4f}px'.rstrip('0').rstrip('.') + 'px')
        except:
            pass
    
    # Scale font-size in style
    if e.get('style') is not None:
        style_str = e.get('style')
        new_style = re.sub(r'font-size:\s*([\d.]+)px', lambda m: f'font-size: {float(m.group(1)) * scale:.2f}px', style_str)
        e.set('style', new_style)
    
    return e

# Add each column
for gi, group in enumerate(groups):
    group_ys = [y for y, elem, orig, short in group]
    min_y = min(group_ys) - 15
    max_y = max(group_ys) + 15
    y_offset = min_y - 30
    
    col_x = MARGIN_L + gi * (col_w + col_gap)
    
    # Grid lines
    for y_val, elem in grid_polylines:
        if min_y <= y_val <= max_y:
            new_elem = shift_scale_shorten(elem, y_offset, col_x, base_scale)
            new_root.append(new_elem)
    
    # Heatmap rects
    for y_val, elem in heatmap_rects:
        if min_y <= y_val <= max_y:
            new_elem = shift_scale_shorten(elem, y_offset, col_x, base_scale)
            new_root.append(new_elem)
    
    # OTU labels (with shortened text)
    for y_val, elem, orig, short in group:
        new_elem = shift_scale_shorten(elem, y_offset, col_x, base_scale, new_text=short)
        new_root.append(new_elem)

# Add column names at bottom center of page
col_name_y = A4_H - MARGIN_B - 5
for elem in col_labels:
    new_elem = copy.deepcopy(elem)
    trans = new_elem.get('transform', '')
    match = re.search(r'translate\(([\d.]+),\s*([\d.]+)\)', trans)
    if match:
        orig_tx = float(match.group(1))
        # Map original x (223, 258, 292, 326) to center of page
        # Original range ~223-326, span ~103, centered at ~274.5
        # Map to page center: A4_W/2 = 421
        # Original relative to center: orig_tx - 274.5
        # New position: 421 + (orig_tx - 274.5) * base_scale
        new_tx = A4_W / 2 + (orig_tx - 274.5) * base_scale
        new_ty = col_name_y
        new_trans = re.sub(r'translate\(([\d.]+),\s*([\d.]+)\)',
                          f'translate({new_tx:.2f},{new_ty:.2f})',
                          trans)
        new_elem.set('transform', new_trans)
        if new_elem.get('style'):
            new_style = re.sub(r'font-size:\s*([\d.]+)px',
                              lambda m: f'font-size: {float(m.group(1)) * base_scale:.2f}px',
                              new_elem.get('style'))
            new_elem.set('style', new_style)
        new_root.append(new_elem)

# Add colorbar top-right
cbar_x = A4_W - MARGIN_R - 70
cbar_y = MARGIN_T + 10

if title_text is not None:
    new_title = copy.deepcopy(title_text)
    new_title.set('x', str(cbar_x))
    new_title.set('y', str(cbar_y))
    if new_title.get('style'):
        new_style = re.sub(r'font-size:\s*([\d.]+)px',
                          lambda m: f'font-size: {float(m.group(1)) * base_scale:.2f}px',
                          new_title.get('style'))
        new_title.set('style', new_style)
    new_root.append(new_title)

if cbar_image is not None:
    new_img = copy.deepcopy(cbar_image)
    img_orig_y = float(cbar_image.get('y', 876))
    img_orig_x = float(cbar_image.get('x', 358))
    x_shift = cbar_x - 358 * base_scale
    y_shift = (cbar_y + 6 + (876 - 869) * base_scale) - 876 * base_scale
    new_img.set('x', str(img_orig_x * base_scale + x_shift))
    new_img.set('y', str(img_orig_y * base_scale + y_shift))
    new_img.set('width', str(float(cbar_image.get('width')) * base_scale))
    new_img.set('height', str(float(cbar_image.get('height')) * base_scale))
    new_root.append(new_img)

for elem in cbar_polylines:
    new_elem = copy.deepcopy(elem)
    pts = new_elem.get('points', '')
    coords = re.findall(r'([\d.]+),([\d.]+)', pts)
    if coords:
        x_shift = cbar_x - 358 * base_scale
        y_shift = (cbar_y + 6 + (876 - 869) * base_scale) - 876 * base_scale
        new_pts = []
        for x_str, y_str in coords:
            new_x = float(x_str) * base_scale + x_shift
            new_y = float(y_str) * base_scale + y_shift
            new_pts.append(f"{new_x:.4f},{new_y:.4f}".rstrip('0').rstrip('.').replace(',.', ',0.').replace('.,', '.0,'))
        new_elem.set('points', ' '.join(new_pts))
        new_root.append(new_elem)

for content, elem in cbar_texts:
    new_elem = copy.deepcopy(elem)
    orig_x = float(new_elem.get('x', 380))
    orig_y = float(new_elem.get('y', 0))
    x_shift = cbar_x - 358 * base_scale
    y_shift = (cbar_y + 6 + (876 - 869) * base_scale) - 876 * base_scale
    new_elem.set('x', str(orig_x * base_scale + x_shift))
    new_elem.set('y', str(orig_y * base_scale + y_shift))
    if new_elem.get('style'):
        new_style = re.sub(r'font-size:\s*([\d.]+)px',
                          lambda m: f'font-size: {float(m.group(1)) * base_scale:.2f}px',
                          new_elem.get('style'))
        new_elem.set('style', new_style)
    new_root.append(new_elem)

# Write output
out_file = "H20_feature_importance_heatmap_plot属_A4_landscape_single_page.svg"
new_tree = ET.ElementTree(new_root)
new_tree.write(out_file, encoding='UTF-8', xml_declaration=True)
print(f"Wrote {out_file}")
print("Done!")
