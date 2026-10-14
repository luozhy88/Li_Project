import xml.etree.ElementTree as ET
import re
import numpy as np

tree = ET.parse('circular_bar_Genus_top50_grouped_from_csv.svg')
root = tree.getroot()
ns = {'svg': 'http://www.w3.org/2000/svg'}

# Find all text elements and locate Lactococcus
texts = root.findall('.//svg:text', ns)
lacto_pos = None
for t in texts:
    for child in t.iter():
        if child.text and 'Lactococcus' in child.text:
            # Get transform
            transform = t.get('transform', '')
            m = re.search(r'translate\(([^)]+)\)', transform)
            if m:
                coords = m.group(1).replace(',', ' ').split()
                x, y = float(coords[0]), float(coords[1])
                lacto_pos = (x, y)
                print(f"Lactococcus label at: ({x:.2f}, {y:.2f})")
            break
    if lacto_pos:
        break

# The chart center should be around (figsize*DPI/2) but in SVG it's in points
# Let's find all paths and check their distance from center
paths = root.findall('.//svg:path', ns)
center_x, center_y = 720, 720  # approximate for 20x20 inch figure at default SVG scaling

# Actually let's compute center from the polar axis circle
# Find paths that look like bars (filled paths with specific colors)
nearby_paths = []
for p in paths:
    d = p.get('d', '')
    style = p.get('style', '')
    fill_match = re.search(r'fill:\s*#([0-9A-Fa-f]{6})', style)
    if not fill_match:
        continue
    fill = fill_match.group(1)
    if fill not in ['e74c3c', '3498db', '2ecc71', '9b59b6']:
        continue
    
    # Parse path for points
    nums = re.findall(r'[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?', d)
    if len(nums) < 4:
        continue
    
    # Take the first point as a proxy for angle/direction
    px, py = float(nums[0]), float(nums[1])
    dist = np.sqrt((px - lacto_pos[0])**2 + (py - lacto_pos[1])**2)
    
    # Compute distance from approximate center to get radius
    center_approx = (720, 720)
    r = np.sqrt((px - center_approx[0])**2 + (py - center_approx[1])**2)
    
    if dist < 200:  # near Lactococcus label
        nearby_paths.append({
            'fill': fill,
            'x': px, 'y': py,
            'dist_to_label': dist,
            'radius': r
        })

print(f"\nFound {len(nearby_paths)} colored paths near Lactococcus")
for p in nearby_paths:
    print(f"  fill=#{p['fill']}, radius={p['radius']:.1f}")
