#!/usr/bin/env python3
"""
Build complete GHRH(1-29) structure.
Strategy:
1. Use PyMOL fab to build the full sequence in extended conformation
2. Align N-terminal residues (1-28) to the 7CZ5 extracted fragment
3. Save the aligned structure
"""

import os, sys
pymol_path = "/home/scroll/miniforge3/lib/python3.12/site-packages"
if os.path.exists(pymol_path) and pymol_path not in sys.path:
    sys.path.insert(0, pymol_path)

import pymol
from pymol import cmd
pymol.finish_launching(["pymol", "-cq"])

BASE = "/home/scroll/personal/cjc-1295/workspace"
STEP1 = f"{BASE}/step1"

# GHRH(1-29) sequence
seq = "YADAIFTNSYRKVLGQLSARKLLQDIMSR"
print(f"Building GHRH(1-29): {seq}")

cmd.reinitialize()
cmd.fab(seq, "GHRH_build", ss=0)

# Load the 7CZ5 extracted fragment for alignment
cmd.load(f"{STEP1}/GHRH_1-29_from_7CZ5.pdb", "GHRH_7CZ5")

# Align built structure to 7CZ5 fragment (first 28 residues)
# Note: 7CZ5 fragment only has 28 residues, so we align resi 1-28
cmd.align("GHRH_build and resi 1-28 and name CA+C+N+O", 
          "GHRH_7CZ5 and name CA+C+N+O")

# Save the aligned full-length structure
cmd.save(f"{STEP1}/GHRH_1-29.pdb", "GHRH_build")

# Verify sequence
cmd.iterate_state(1, "GHRH_build and name CA", "print('  ',resi,resn)")

print(f"\nSaved: {STEP1}/GHRH_1-29.pdb")

# Also create a version with N-terminus NH3+ and C-terminus COO- explicitly
# PyMOL fab does this by default, but let's make sure termini are correct

cmd.quit()
