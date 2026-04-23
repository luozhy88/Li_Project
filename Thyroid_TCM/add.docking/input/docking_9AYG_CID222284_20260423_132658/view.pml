# Auto-generated PyMOL script
load /tmp/RtmpL6MDX0/docking_20260423_132439/receptor/receptor_clean.pdb
load /tmp/RtmpL6MDX0/docking_20260423_132439/best_pose.sdf

remove solvent
select ligand, organic
select protein, all and not ligand

hide everything
bg_color white
show cartoon, protein
color gray80, protein
show sticks, ligand
util.cbag ligand

select contacts, byres protein within 5 of ligand
show sticks, contacts
util.cbag contacts

distance hbonds, ligand, protein, mode=2, cutoff=3.5
hide labels, hbonds
color yellow, hbonds
set dash_gap, 0.3
set dash_radius, 0.05

label contacts and n. CA, "%s-%s" % (resn, resi)
set label_size, 10
set label_font_id, 7
set label_position, (0,0,2)
set ray_opaque_background, 0

orient ligand
zoom ligand, 5
png /tmp/RtmpL6MDX0/docking_20260423_132439/output_local.png, width=2400, dpi=300, ray=1
save /tmp/RtmpL6MDX0/docking_20260423_132439/output_local.pov

orient
png /tmp/RtmpL6MDX0/docking_20260423_132439/output_global.png, width=2400, dpi=300, ray=1
save /tmp/RtmpL6MDX0/docking_20260423_132439/output_global.pov

save /tmp/RtmpL6MDX0/docking_20260423_132439/session.pse
quit
