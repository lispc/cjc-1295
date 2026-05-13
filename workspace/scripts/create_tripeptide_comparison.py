#!/usr/bin/env python3
"""
Create a comparison figure of WT YAD vs D-Ala2 YDD in the DPP-IV active site.
Uses PyMOL via pml script generation.
"""

import os
import subprocess

BASE = "/home/scroll/personal/cjc-1295/workspace"
STEP0 = f"{BASE}/step0"
STEP1 = f"{BASE}/step1"
RESULTS = f"{BASE}/results"

os.makedirs(RESULTS, exist_ok=True)

# Write PyMOL script
pml_script = f"""
# Load structures
load {STEP1}/DPP4_clean.pdb, DPP4
load {STEP0}/YAD_aligned_to_DPP4.pdb, YAD_WT
load {STEP0}/YDD_tripeptide.pdb, YAD_DA

# Align D-Ala2 YDD to WT YAD
pair_fit YAD_DA and name N+CA+C+O, YAD_WT and name N+CA+C+O

# Create complex for visualization
create complex_WT, DPP4 or YAD_WT
create complex_DA, DPP4 or YAD_DA

# Hide everything initially
hide everything, all

# Show DPP-IV as surface with transparency
show surface, DPP4
color white, DPP4
set transparency, 0.6, DPP4

# Show catalytic residues as sticks
show sticks, DPP4 and resi 630+205+206+547

# Color catalytic residues
color red, DPP4 and resi 630 and name OG     # Ser630 nucleophile
color blue, DPP4 and resi 205+206 and name OE1+OE2  # Glu205/206
color orange, DPP4 and resi 547 and name OH   # Tyr547 oxyanion hole

# Show peptides as sticks
show sticks, YAD_WT
show sticks, YAD_DA

# Color WT Ala2
color yellow, YAD_WT and resi 2 and name CA+CB+HB1+HB2+HB3
color green, YAD_WT and resi 2 and name N+C+O+CA

# Color D-Ala2 differently
color cyan, YAD_DA and resi 2 and name CA+CB+HB1+HB2+HB3
color magenta, YAD_DA and resi 2 and name N+C+O+CA

# Color rest of peptides
color gray, YAD_WT and resi 1+3
color gray, YAD_DA and resi 1+3

# Show N-terminus and C-terminus
color blue, YAD_WT and resi 1 and name N
color blue, YAD_DA and resi 1 and name N
color red, YAD_WT and resi 3 and name O
color red, YAD_DA and resi 3 and name O

# Label key residues
label DPP4 and resi 630 and name OG, "Ser630 Oγ"
label DPP4 and resi 205 and name CD, "Glu205"
label DPP4 and resi 206 and name CD, "Glu206"
label DPP4 and resi 547 and name CZ, "Tyr547"
label YAD_WT and resi 2 and name CA, "L-Ala2"
label YAD_DA and resi 2 and name CA, "D-Ala2"

# Set view to active site
center DPP4 and resi 630+205+206+547
zoom DPP4 and resi 630+205+206+547, 8

# Set rendering quality
set ray_trace_mode, 1
set ray_trace_fog, 0
set antialias, 2
set ambient, 0.4
set reflect, 0.5
set specular, 0.5

# First view: WT only
hide sticks, YAD_DA
hide labels, YAD_DA
show labels, YAD_WT
show labels, DPP4 and resi 630+205+206+547
ray 2400, 2400
png {RESULTS}/step0_WT_YAD_active_site.png, dpi=300, width=2400, height=2400

# Second view: D-Ala2 only
hide sticks, YAD_WT
hide labels, YAD_WT
show sticks, YAD_DA
show labels, YAD_DA
show labels, DPP4 and resi 630+205+206+547
ray 2400, 2400
png {RESULTS}/step0_DA_YAD_active_site.png, dpi=300, width=2400, height=2400

# Third view: Overlay comparison
show sticks, YAD_WT
show labels, YAD_WT
show labels, YAD_DA
set label_position, (0,0,0)

# Different angle for overlay
rotate y, 15
rotate x, 10
ray 2400, 2400
png {RESULTS}/step0_overlay_comparison.png, dpi=300, width=2400, height=2400

# Save session
save {RESULTS}/step0_comparison.pse

quit
"""

pml_path = f"{BASE}/scripts/tripeptide_comparison.pml"
with open(pml_path, "w") as f:
    f.write(pml_script)

print(f"PyMOL script written to {pml_path}")

# Run PyMOL
result = subprocess.run(
    ["/home/scroll/miniforge3/bin/pymol", "-cq", pml_path],
    capture_output=True, text=True
)

if result.returncode == 0:
    print("PyMOL rendering completed successfully")
    print(f"  {RESULTS}/step0_WT_YAD_active_site.png")
    print(f"  {RESULTS}/step0_DA_YAD_active_site.png")
    print(f"  {RESULTS}/step0_overlay_comparison.png")
else:
    print(f"PyMOL failed: {result.stderr}")
    print(f"stdout: {result.stdout[-500:]}")
