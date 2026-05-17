import xml.etree.ElementTree as ET
import copy
import re

# Parse SVG
svg_file = "H20_feature_importance_heatmap_plot属.svg"
tree = ET.parse(svg_file)
root = tree.getroot()

ns = {'svg': 'http://www.w3.org/2000/svg'}
ET.register_namespace('', 'http://www.w3.org/2000/svg')
ET.register_namespace('xlink', 'http://www.w3.org/1999/xlink')

# Extract OTU label text elements
otu_texts = []
for elem in root.iter():
    if elem.tag == '{http://www.w3.org/2000/svg}text' and elem.get('text-anchor') == 'end' and elem.get('x') == '194.15':
        y = float(elem.get('y', 0))
        otu_texts.append((y, elem))

# Sort by y ascending (top to bottom)
otu_texts.sort(key=lambda x: x[0])
print(f"Total OTU rows: {len(otu_texts)}")

# Group into 4 segments
n_groups = 4
rows_per_group = len(otu_texts) // n_groups
remainder = len(otu_texts) % n_groups

groups = []
start = 0
for i in range(n_groups):
    size = rows_per_group + (1 if i < remainder else 0)
    group_texts = otu_texts[start:start+size]
    groups.append(group_texts)
    start += size

print(f"Groups: {[len(g) for g in groups]}")

# Collect all rects by y
all_rects = []
for elem in root.iter():
    if elem.tag == '{http://www.w3.org/2000/svg}rect':
        try:
            y = float(elem.get('y', 0))
            width = elem.get('width', '')
            all_rects.append((y, width, elem))
        except:
            pass

heatmap_rects = [(y, e) for y, w, e in all_rects if w == '33.92']

# Collect all polylines by y
all_polylines = []
for elem in root.iter():
    if elem.tag == '{http://www.w3.org/2000/svg}polyline':
        pts = elem.get('points', '')
        coords = re.findall(r'[\d.]+', pts)
        if len(coords) >= 2:
            y_val = float(coords[1])
            all_polylines.append((y_val, elem))

# Collect column labels
col_labels = []
for elem in root.iter():
    if elem.tag == '{http://www.w3.org/2000/svg}text' and 'transform' in elem.attrib and 'rotate' in elem.attrib.get('transform', ''):
        col_labels.append(elem)

# Collect colorbar elements
cbar_image = None
cbar_polylines = []
cbar_texts = []  # 0.00, 0.25, etc.
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

# Sort cbar texts by content order
order = {'0.00': 0, '0.25': 1, '0.50': 2, '0.75': 3, '1.00': 4}
cbar_texts.sort(key=lambda x: order.get(x[0], 99))

# For each group, build a new SVG
for gi, group in enumerate(groups):
    group_ys = [y for y, t in group]
    min_y = min(group_ys) - 15  # padding top
    max_y = max(group_ys) + 15  # padding bottom
    
    # Create new root
    new_root = ET.Element('svg')
    new_root.set('class', 'svglite')
    
    # Calculate new height
    content_height = max_y - min_y
    new_height = content_height + 40  # extra margin for column labels
    new_width = 432.00
    
    new_root.set('width', f'{new_width:.2f}pt')
    new_root.set('height', f'{new_height:.2f}pt')
    new_root.set('viewBox', f'0 0 {new_width:.2f} {new_height:.2f}')
    
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
    
    # Add background rect
    bg = ET.SubElement(new_root, 'rect')
    bg.set('width', '100%')
    bg.set('height', '100%')
    bg.set('style', 'stroke: none; fill: #FFFFFF;')
    
    # y offset: shift content so that min_y maps to 30pt from top
    y_offset = min_y - 30
    
    # Helper to shift element y
    def shift_element_y(elem, y_offset):
        e = copy.deepcopy(elem)
        if e.get('y') is not None:
            try:
                old_y = float(e.get('y'))
                e.set('y', f'{old_y - y_offset:.6f}'.rstrip('0').rstrip('.'))
            except:
                pass
        if e.get('points') is not None:
            pts = e.get('points')
            coords = re.findall(r'([\d.]+),([\d.]+)', pts)
            new_pts = []
            for x_str, y_str in coords:
                new_y = float(y_str) - y_offset
                new_pts.append(f"{x_str},{new_y:.6f}".rstrip('0').rstrip('.'))
            e.set('points', ' '.join(new_pts))
        if e.get('transform') is not None:
            trans = e.get('transform')
            new_trans = re.sub(r'translate\(([\d.]+),\s*([\d.]+)\)', 
                              lambda m: f'translate({m.group(1)},{float(m.group(2)) - y_offset:.2f})',
                              trans)
            e.set('transform', new_trans)
        return e
    
    # Helper to place element at fixed y
    def fix_element_y(elem, target_y):
        e = copy.deepcopy(elem)
        if e.get('y') is not None:
            try:
                old_y = float(e.get('y'))
                e.set('y', f'{target_y:.6f}'.rstrip('0').rstrip('.'))
            except:
                pass
        if e.get('points') is not None:
            pts = e.get('points')
            coords = re.findall(r'([\d.]+),([\d.]+)', pts)
            new_pts = []
            for x_str, y_str in coords:
                # shift by same amount as first point
                first_y = float(coords[0][1]) if coords else 0
                dy = target_y - first_y
                new_y = float(y_str) + dy
                new_pts.append(f"{x_str},{new_y:.6f}".rstrip('0').rstrip('.'))
            e.set('points', ' '.join(new_pts))
        if e.get('transform') is not None:
            trans = e.get('transform')
            match = re.search(r'translate\(([\d.]+),\s*([\d.]+)\)', trans)
            if match:
                tx, ty = float(match.group(1)), float(match.group(2))
                new_ty = target_y  # place at target_y
                new_trans = re.sub(r'translate\(([\d.]+),\s*([\d.]+)\)', 
                                  f'translate({tx:.2f},{new_ty:.2f})',
                                  trans)
                e.set('transform', new_trans)
        return e
    
    # Add heatmap polylines (grid lines) within y range
    for y_val, elem in all_polylines:
        if min_y <= y_val <= max_y and elem not in cbar_polylines:
            new_elem = shift_element_y(elem, y_offset)
            new_root.append(new_elem)
    
    # Add heatmap rects within y range
    for y_val, elem in heatmap_rects:
        if min_y <= y_val <= max_y:
            new_elem = shift_element_y(elem, y_offset)
            new_root.append(new_elem)
    
    # Add OTU texts
    for y_val, elem in group:
        new_elem = shift_element_y(elem, y_offset)
        new_root.append(new_elem)
    
    # Add column labels at fixed bottom position (last group only, near bottom of its view)
    if gi == n_groups - 1:
        col_label_y = new_height - 20  # 20pt from bottom
        for elem in col_labels:
            new_elem = fix_element_y(elem, col_label_y)
            new_root.append(new_elem)
    
    # Add colorbar, title, and scale labels at fixed top-right position (first group only)
    if gi == 0:
        cbar_top_y = 50  # 50pt from top
        # Original relative positions:
        # title at y=869.38, image at y=876, ticks from 876 to 962, texts from 879 to 965
        # We set title y to cbar_top_y
        if title_text is not None:
            new_elem = copy.deepcopy(title_text)
            new_elem.set('y', str(cbar_top_y))
            new_root.append(new_elem)
        
        # Image: original y=876, title y=869.38, diff = 6.62
        if cbar_image is not None:
            new_elem = copy.deepcopy(cbar_image)
            new_elem.set('y', str(cbar_top_y + 6.62))
            new_root.append(new_elem)
        
        # Cbar polylines: original y from 876.15 to 962.26
        for elem in cbar_polylines:
            new_elem = copy.deepcopy(elem)
            pts = new_elem.get('points', '')
            coords = re.findall(r'([\d.]+),([\d.]+)', pts)
            if coords:
                first_y = float(coords[0][1])
                dy = (cbar_top_y + 6.62) - 876.00  # approximate original image y
                new_pts = []
                for x_str, y_str in coords:
                    new_y = float(y_str) + dy
                    new_pts.append(f"{x_str},{new_y:.6f}".rstrip('0').rstrip('.'))
                new_elem.set('points', ' '.join(new_pts))
            new_root.append(new_elem)
        
        # Cbar texts: original y from 879.36 to 965.47
        for content, elem in cbar_texts:
            new_elem = copy.deepcopy(elem)
            try:
                old_y = float(new_elem.get('y', 0))
                dy = (cbar_top_y + 6.62) - 876.00
                new_y = old_y + dy
                new_elem.set('y', f'{new_y:.6f}'.rstrip('0').rstrip('.'))
            except:
                pass
            new_root.append(new_elem)
    
    # Write output
    out_file = f"H20_feature_importance_heatmap_plot属_part{gi+1}of{n_groups}.svg"
    new_tree = ET.ElementTree(new_root)
    new_tree.write(out_file, encoding='UTF-8', xml_declaration=True)
    print(f"Wrote {out_file}: {new_height:.0f}pt height")

print("Done!")
