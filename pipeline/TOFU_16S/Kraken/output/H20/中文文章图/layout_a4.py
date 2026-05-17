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
        otu_texts.append((y, elem))

otu_texts.sort(key=lambda x: x[0])
n_rows = len(otu_texts)
print(f"Total OTU rows: {n_rows}")

# Group into 4 panels (top to bottom order in original)
n_groups = 4
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
            # Exclude colorbar polylines
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

# Helper functions
def get_y_range(elem):
    """Get y coordinate(s) from element"""
    if elem.get('y') is not None:
        try:
            return [float(elem.get('y'))]
        except:
            pass
    if elem.get('points') is not None:
        pts = elem.get('points')
        coords = re.findall(r'[\d.]+', pts)
        if len(coords) >= 2:
            return [float(coords[i]) for i in range(1, len(coords), 2)]
    if elem.get('transform') is not None:
        trans = elem.get('transform')
        match = re.search(r'translate\(([\d.]+),\s*([\d.]+)\)', trans)
        if match:
            return [float(match.group(2))]
    return []

def shift_and_scale_element(elem, y_offset, x_offset, scale):
    """Apply y shift, x shift, and uniform scale to element"""
    e = copy.deepcopy(elem)
    
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
            new_x = (old_x) * scale + x_offset
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
    
    # Scale textLength
    if e.get('textLength') is not None:
        try:
            tl = e.get('textLength').replace('px', '')
            new_tl = float(tl) * scale
            e.set('textLength', f'{new_tl:.4f}px'.rstrip('0').rstrip('.') + 'px')
        except:
            pass
    
    # Scale font-size in style
    if e.get('style') is not None:
        style = e.get('style')
        new_style = re.sub(r'font-size:\s*([\d.]+)px', lambda m: f'font-size: {float(m.group(1)) * scale:.2f}px', style)
        e.set('style', new_style)
    
    return e

# A4 Landscape layout parameters
A4_W, A4_H = 842, 595
MARGIN_L, MARGIN_R = 50, 50
MARGIN_T, MARGIN_B = 50, 50
usable_w = A4_W - MARGIN_L - MARGIN_R
usable_h = A4_H - MARGIN_T - MARGIN_B

panel_gap = 40
panel_w = (usable_w - panel_gap) / 2

# Determine scale from original panel width (approx 432pt)
# But we want to fit within panel_w while preserving aspect
# Original content width: text starts near x=0, heatmap ends at x~342, colorbar at x~380
# Let's use full original width 432 as reference
base_scale = panel_w / 432.0
print(f"Panel scale: {base_scale:.3f}")

# Original row height ~9.1pt, 49 rows -> ~446pt content height
# Scaled height = 446 * base_scale
content_h = 500 * base_scale  # approximate panel content height
panel_v_margin = (usable_h - content_h) / 2
print(f"Panel content height: {content_h:.1f}pt, vertical margin: {panel_v_margin:.1f}pt")

# Build 2 pages
pages = [
    (0, 1),  # page 1: groups 0 and 1
    (2, 3),  # page 2: groups 2 and 3
]

for page_idx, (g1, g2) in enumerate(pages):
    # Create A4 landscape SVG
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
    
    # Helper to add a group panel
    def add_panel(gi, panel_x_offset):
        group = groups[gi]
        group_ys = [y for y, t in group]
        min_y = min(group_ys) - 15
        max_y = max(group_ys) + 15
        
        y_offset = min_y - 30
        x_offset = MARGIN_L + panel_x_offset
        
        # Grid lines
        for y_val, elem in grid_polylines:
            if min_y <= y_val <= max_y:
                new_elem = shift_and_scale_element(elem, y_offset, x_offset, base_scale)
                new_root.append(new_elem)
        
        # Heatmap rects
        for y_val, elem in heatmap_rects:
            if min_y <= y_val <= max_y:
                new_elem = shift_and_scale_element(elem, y_offset, x_offset, base_scale)
                new_root.append(new_elem)
        
        # OTU labels
        for y_val, elem in group:
            new_elem = shift_and_scale_element(elem, y_offset, x_offset, base_scale)
            new_root.append(new_elem)
        
        # Column labels at bottom of panel
        col_label_y = MARGIN_T + usable_h - 30  # base y for bottom area
        # Original column labels are at transform y ~1808. We need to place them near panel bottom
        for elem in col_labels:
            new_elem = copy.deepcopy(elem)
            trans = new_elem.get('transform', '')
            match = re.search(r'translate\(([\d.]+),\s*([\d.]+)\)', trans)
            if match:
                orig_tx = float(match.group(1))
                # Scale x relative to heatmap center
                new_tx = orig_tx * base_scale + x_offset
                new_ty = col_label_y
                new_trans = re.sub(r'translate\(([\d.]+),\s*([\d.]+)\)',
                                  f'translate({new_tx:.2f},{new_ty:.2f})',
                                  trans)
                new_elem.set('transform', new_trans)
                # Scale font-size
                if new_elem.get('style'):
                    new_style = re.sub(r'font-size:\s*([\d.]+)px',
                                      lambda m: f'font-size: {float(m.group(1)) * base_scale:.2f}px',
                                      new_elem.get('style'))
                    new_elem.set('style', new_style)
                new_root.append(new_elem)
    
    # Add left panel
    add_panel(g1, 0)
    # Add right panel
    add_panel(g2, panel_w + panel_gap)
    
    # Add colorbar on page 1 only, top-right area but clear of right panel
    if page_idx == 0 and cbar_image is not None:
        # Place colorbar at x=750 (right of right panel heatmap which ends at ~718)
        cbar_x = A4_W - MARGIN_R - 40
        cbar_y = MARGIN_T + 10
        
        # Title
        if title_text is not None:
            new_title = copy.deepcopy(title_text)
            new_title.set('x', str(cbar_x))
            new_title.set('y', str(cbar_y))
            # Scale font
            if new_title.get('style'):
                new_style = re.sub(r'font-size:\s*([\d.]+)px',
                                  lambda m: f'font-size: {float(m.group(1)) * base_scale:.2f}px',
                                  new_title.get('style'))
                new_title.set('style', new_style)
            new_root.append(new_title)
        
        # Image
        new_img = copy.deepcopy(cbar_image)
        img_orig_y = float(cbar_image.get('y', 876))
        img_orig_x = float(cbar_image.get('x', 358))
        new_img.set('x', str(cbar_x + (img_orig_x - 358) * base_scale))
        new_img.set('y', str(cbar_y + 6 + (img_orig_y - 869) * base_scale))
        new_img.set('width', str(float(cbar_image.get('width')) * base_scale))
        new_img.set('height', str(float(cbar_image.get('height')) * base_scale))
        new_root.append(new_img)
        
        # Tick polylines
        for elem in cbar_polylines:
            new_elem = copy.deepcopy(elem)
            pts = new_elem.get('points', '')
            coords = re.findall(r'([\d.]+),([\d.]+)', pts)
            if coords:
                first_y = float(coords[0][1])
                dy = (cbar_y + 6 + (876 - 869) * base_scale) - 876 * base_scale
                # Actually simpler: compute offset from original image position
                orig_img_y = 876
                new_img_y = cbar_y + 6 + (876 - 869) * base_scale
                y_shift = new_img_y - orig_img_y * base_scale
                x_shift = cbar_x - 358 * base_scale
                new_pts = []
                for x_str, y_str in coords:
                    new_x = float(x_str) * base_scale + x_shift
                    new_y = float(y_str) * base_scale + y_shift
                    new_pts.append(f"{new_x:.4f},{new_y:.4f}".rstrip('0').rstrip('.').replace(',.', ',0.').replace('.,', '.0,'))
                new_elem.set('points', ' '.join(new_pts))
                new_root.append(new_elem)
        
        # Tick labels
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
    out_file = f"H20_feature_importance_heatmap_plot属_A4_landscape_page{page_idx+1}.svg"
    new_tree = ET.ElementTree(new_root)
    new_tree.write(out_file, encoding='UTF-8', xml_declaration=True)
    print(f"Wrote {out_file}")

print("Done!")
