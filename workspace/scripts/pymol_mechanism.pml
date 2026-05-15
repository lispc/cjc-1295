# PyMOL script: D-Ala2 vs L-Ala2 catalytic geometry at DPP-IV active site
# Shows the "lock but wrong angle" inhibition mechanism

# ── Setup ──
bg_color white
set ray_trace_mode, 1
set antialias, 2
set cartoon_fancy_helices, 1
set depth_cue, 0

# ── Load structures ──
load frame_DALA_100ns.pdb, DALA_complex
load frame_LALA_25ns.pdb, LALA_complex
load frame_start.pdb, start_complex

# ── Color setup ──
set_color dppiv_grey  = [0.75, 0.75, 0.78]
set_color catalytic   = [0.90, 0.45, 0.20]   # orange
set_color oxyanion    = [0.18, 0.65, 0.55]    # teal
set_color d_ala_color = [0.90, 0.25, 0.30]    # red
set_color l_ala_color = [0.28, 0.47, 0.62]    # blue

# ────────────────────────────────────────────────────────
# LEFT PANEL: D-Ala2 (CJC-1295) — tight but wrong angle
# ────────────────────────────────────────────────────────
viewport 800, 800
set viewport 0, 400, 0, 800

# Load D-Ala complex
disable all
enable DALA_complex

# DPP-IV: light cartoon
show cartoon, DALA_complex and chain A
color dppiv_grey, DALA_complex and chain A
set cartoon_transparency, 0.65, DALA_complex and chain A

# Catalytic triad: sticks
select catalytic_res, (DALA_complex and chain A and resi 630+708+740)
show sticks, catalytic_res
color catalytic, catalytic_res
set stick_radius, 0.2, catalytic_res

# Oxyanion hole
select oxy_res, (DALA_complex and chain A and resi 547+631)
show sticks, oxy_res
color oxyanion, oxy_res

# Glu205/206 anchor
select anchor_res, (DALA_complex and chain A and resi 205+206)
show sticks, anchor_res
color tv_red, anchor_res

# GHRH peptide (chain B) — thick sticks
show sticks, DALA_complex and chain B and resi 1-5
color d_ala_color, DALA_complex and chain B and resi 1-5
set stick_radius, 0.3, DALA_complex and chain B

# Highlight Ala2 sidechain
show spheres, DALA_complex and chain B and resi 2 and name CB
color d_ala_color, DALA_complex and chain B and resi 2 and name CB
set sphere_scale, 0.5, DALA_complex and chain B and resi 2 and name CB

# Distance Ser630 OG → Ala2 C
distance dala_dist, (DALA_complex and chain A and resi 630 and name OG), (DALA_complex and chain B and resi 2 and name C)
color d_ala_color, dala_dist
set dash_length, 0.15, dala_dist
set dash_gap, 0.1, dala_dist

# Angle Ser630 OG – Ala2 C – Ala2 N
angle dala_angle, (DALA_complex and chain A and resi 630 and name OG), (DALA_complex and chain B and resi 2 and name C), (DALA_complex and chain B and resi 2 and name N)
color d_ala_color, dala_angle

# Labels
set label_size, 20
set label_color, d_ala_color
set label_font_id, 5

# Zoom on active site
zoom (DALA_complex and chain A and resi 630), 12
clip near, 0
clip far, 40

# Scene label
pseudoatom dala_label, pos=[0,0,0]
set_label dala_label, "D-Ala2"

# ────────────────────────────────────────────────────────
# RIGHT PANEL: L-Ala2 (native GHRH) — drifted away
# ────────────────────────────────────────────────────────
set viewport 400, 800, 0, 800

disable all
enable LALA_complex

# DPP-IV: light cartoon
show cartoon, LALA_complex and chain A
color dppiv_grey, LALA_complex and chain A
set cartoon_transparency, 0.65, LALA_complex and chain A

# Catalytic triad
select catalytic_res2, (LALA_complex and chain A and resi 630+708+740)
show sticks, catalytic_res2
color catalytic, catalytic_res2

# Oxyanion hole
select oxy_res2, (LALA_complex and chain A and resi 547+631)
show sticks, oxy_res2
color oxyanion, oxy_res2

# Glu205/206
select anchor_res2, (LALA_complex and chain A and resi 205+206)
show sticks, anchor_res2
color tv_red, anchor_res2

# GHRH peptide
show sticks, LALA_complex and chain B and resi 1-5
color l_ala_color, LALA_complex and chain B and resi 1-5
set stick_radius, 0.3, LALA_complex and chain B

# Highlight Ala2
show spheres, LALA_complex and chain B and resi 2 and name CB
color l_ala_color, LALA_complex and chain B and resi 2 and name CB
set sphere_scale, 0.5, LALA_complex and chain B and resi 2 and name CB

# Distance
distance lala_dist, (LALA_complex and chain A and resi 630 and name OG), (LALA_complex and chain B and resi 2 and name C)
color l_ala_color, lala_dist

# Angle
angle lala_angle, (LALA_complex and chain A and resi 630 and name OG), (LALA_complex and chain B and resi 2 and name C), (LALA_complex and chain B and resi 2 and name N)
color l_ala_color, lala_angle

# Same zoom
zoom (LALA_complex and chain A and resi 630), 12

# ── Ray-trace and save ──
set ray_opaque_background, off
viewport 1600, 800
ray 1600, 800
png mechanism_pymol.png, dpi=300

print "Saved mechanism_pymol.png"
