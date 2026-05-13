#!/usr/bin/env python3
"""
Re-align YAD tripeptide to Diprotin A using pair_fit.
Assign chain Y to YAD to avoid residue name collisions with DPP-IV.
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

cmd.reinitialize()
cmd.load(f"{STEP1}/DPP4_with_diprotinA.pdb", "ref")
cmd.load(f"{STEP0}/YAD_tripeptide.pdb", "yad")

# Assign chain Y to YAD before alignment
cmd.alter("yad", "chain='Y'")
cmd.sort()

# Use pair_fit to superimpose YAD onto Diprotin A (chain D)
cmd.pair_fit(
    "yad and resi 1 and name N+CA+C+O", "ref and chain D and resi 1 and name N+CA+C+O",
    "yad and resi 2 and name N+CA+C+O", "ref and chain D and resi 2 and name N+CA+C+O",
    "yad and resi 3 and name N+CA+C+O", "ref and chain D and resi 3 and name N+CA+C+O"
)

# Save aligned YAD alone
cmd.save(f"{STEP0}/YAD_aligned_to_DPP4.pdb", "yad")

# Save complex with YAD and DPP-IV (remove Diprotin A chain D)
cmd.remove("chain D")
cmd.create("complex", "ref or yad")
cmd.save(f"{STEP0}/DPP4_with_YAD.pdb", "complex")

print("Re-alignment complete with chain Y assigned to YAD.")
print(f"  Saved: {STEP0}/YAD_aligned_to_DPP4.pdb")
print(f"  Saved: {STEP0}/DPP4_with_YAD.pdb")

cmd.quit()
