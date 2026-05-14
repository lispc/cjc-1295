#!/usr/bin/env python3
"""
Phase 1 Structure Preparation Script (CURRENT — use this one)
- Preprocess 1NU8 (DPP-IV) -> DPP4_clean.pdb
  IMPORTANT: Diprotin A (chain D) binds to chain B, not chain A!
  We keep chain B as the receptor.
- Extract GHRH(1-29) from 7CZ5 -> GHRH_1-29.pdb
- Build YAD tripeptide -> YAD_tripeptide.pdb

NOTE: prepare_structures.py (v1) is DEPRECATED. It keeps chain A, which is wrong.
      This v2 is the correct version. Always use v2 for new runs.
"""

import os, sys

pymol_path = "/home/scroll/miniforge3/lib/python3.12/site-packages"
if os.path.exists(pymol_path) and pymol_path not in sys.path:
    sys.path.insert(0, pymol_path)

import pymol
from pymol import cmd
pymol.finish_launching(["pymol", "-cq"])

BASE = "/home/scroll/personal/cjc-1295/workspace"
STRUCT = f"{BASE}/structures"
STEP0 = f"{BASE}/step0"
STEP1 = f"{BASE}/step1"

os.makedirs(STEP0, exist_ok=True)
os.makedirs(STEP1, exist_ok=True)

print("=" * 60)
print("Step 1.1: Preprocessing 1NU8 -> DPP4_clean.pdb")
print("NOTE: Diprotin A (chain D) binds to CHAIN B, keeping chain B")
print("=" * 60)

cmd.reinitialize()
cmd.load(f"{STRUCT}/1NU8.pdb", "1NU8")

# Remove solvent and unwanted heteroatoms
cmd.remove("solvent")
cmd.remove("resn HOH")
cmd.remove("resn SO4")
cmd.remove("resn NAG")

# Keep chain B (the monomer with Diprotin A bound) and chain D (Diprotin A)
# Chain A is the other monomer, remove it
cmd.create("DPP4_with_ligand", "chain B or chain D")
cmd.create("DPP4_clean", "chain B")

# Rename chain B to A for consistency with the rest of the pipeline
cmd.alter("chain B", "chain='A'")
cmd.alter("chain D", "chain='D'")
cmd.sort()

# Save clean receptor
cmd.save(f"{STEP1}/DPP4_clean.pdb", "DPP4_clean")
# Save with ligand for reference
cmd.save(f"{STEP1}/DPP4_with_diprotinA.pdb", "DPP4_with_ligand")

print(f"  Saved: {STEP1}/DPP4_clean.pdb")
print(f"  Saved: {STEP1}/DPP4_with_diprotinA.pdb")

# Verify Diprotin A position relative to catalytic residues
print("\n--- Diprotin A (chain D) position check ---")
cmd.iterate_state(1, "chain D and name CA", "print('  DiprotinA res:',resn,resi)")

# Check distances from Diprotin A to catalytic residues
from chempy import cpv
coords = {}
for sel in ["chain D and resi 2 and name CA", "resi 630 and name CA", "resi 205 and name CA", "resi 206 and name CA"]:
    model = cmd.get_model(sel)
    if model.atom:
        a = model.atom[0]
        coords[sel] = a.coord
        print(f"  {sel}: ({a.coord[0]:.1f}, {a.coord[1]:.1f}, {a.coord[2]:.1f})")

print("\n" + "=" * 60)
print("Step 1.2: Extract GHRH(1-29) from 7CZ5")
print("=" * 60)

cmd.reinitialize()
cmd.load(f"{STRUCT}/7CZ5.pdb", "7CZ5")

# Chain P is GHRH
cmd.create("GHRH_1-29", "chain P and resi 1-29")
cmd.save(f"{STEP1}/GHRH_1-29_from_7CZ5.pdb", "GHRH_1-29")

print("  Extracted GHRH(1-29) from 7CZ5 chain P")
cmd.iterate_state(1, "GHRH_1-29 and name CA", "print('  ',resi,resn)")

print("\n" + "=" * 60)
print("Step 0: Build YAD tripeptide")
print("=" * 60)

cmd.reinitialize()
# CRITICAL FIX: Do NOT use cmd.fab("YAD", ss=0) which builds a linear
# peptide with phi=psi=0 (Ramachandran forbidden). Instead, extract
# the natural YAD conformation from the GHRH crystal structure (7CZ5).
cmd.load(f"{STEP1}/GHRH_1-29_from_7CZ5.pdb", "GHRH_ref")
cmd.create("YAD", "GHRH_ref and resi 1-3")
cmd.delete("GHRH_ref")
cmd.alter("YAD", "chain='Y'")
cmd.sort()
cmd.save(f"{STEP0}/YAD_tripeptide.pdb", "YAD")
print(f"  Saved: {STEP0}/YAD_tripeptide.pdb")
print("  NOTE: YAD extracted from native GHRH(1-29) crystal structure")
print("        (avoids linear fab() artifact with phi=psi=0)")

cmd.quit()
