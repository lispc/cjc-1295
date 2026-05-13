#!/usr/bin/env python3
"""
Prepare the starting structure for Rosetta FlexPepDock.
Align GHRH(1-29) N-terminus to the YAD tripeptide position in DPP-IV.
"""

import os, sys
pymol_path = "/home/scroll/miniforge3/lib/python3.12/site-packages"
if os.path.exists(pymol_path) and pymol_path not in sys.path:
    sys.path.insert(0, pymol_path)

import pymol
from pymol import cmd
pymol.finish_launching(["pymol", "-cq"])

BASE = "/home/scroll/personal/cjc-1295/workspace"
STEP0 = f"{BASE}/step0"
STEP1 = f"{BASE}/step1"
STEP2 = f"{BASE}/step2"

os.makedirs(STEP2, exist_ok=True)

cmd.reinitialize()
cmd.load(f"{STEP1}/DPP4_clean.pdb", "DPP4")
cmd.load(f"{STEP0}/YAD_aligned_to_DPP4.pdb", "YAD_ref")
cmd.load(f"{STEP1}/GHRH_1-29.pdb", "GHRH")

# Assign chain B to GHRH (so receptor is A, peptide is B)
cmd.alter("GHRH", "chain='B'")
cmd.sort()

# Align GHRH residues 1-3 (YAD) to the reference YAD tripeptide
cmd.pair_fit(
    "GHRH and resi 1 and name N+CA+C+O", "YAD_ref and resi 1 and name N+CA+C+O",
    "GHRH and resi 2 and name N+CA+C+O", "YAD_ref and resi 2 and name N+CA+C+O",
    "GHRH and resi 3 and name N+CA+C+O", "YAD_ref and resi 3 and name N+CA+C+O"
)

# Remove YAD reference
cmd.delete("YAD_ref")

# Save the starting complex
cmd.save(f"{STEP2}/DPP4_GHRH_start.pdb", "DPP4 or GHRH")
print(f"Saved: {STEP2}/DPP4_GHRH_start.pdb")

# Verify placement
print("\n--- GHRH N-terminus placement check ---")
cmd.iterate_state(1, "GHRH and resi 1-3 and name CA", "print('  GHRH',resi,resn,chain)")

# Quick distance check
cmd.select("ser630", "DPP4 and resi 630 and name OG")
cmd.select("ala2_c", "GHRH and resi 2 and name C")
cmd.select("glu205", "DPP4 and resi 205 and name OE1")
cmd.select("glu206", "DPP4 and resi 206 and name OE1")

# Use get_distance
d1 = cmd.get_distance("ser630", "ala2_c")
d2 = cmd.get_distance("DPP4 and resi 205 and name CA", "GHRH and resi 1 and name N")
d3 = cmd.get_distance("DPP4 and resi 206 and name CA", "GHRH and resi 1 and name N")

print(f"  Ser630 OG to GHRH Ala2 C: {d1:.2f} Å")
print(f"  Glu205 CA to GHRH Tyr1 N: {d2:.2f} Å")
print(f"  Glu206 CA to GHRH Tyr1 N: {d3:.2f} Å")

cmd.quit()
