#!/usr/bin/env python3
"""
DEPRECATED — Use prepare_structures_v2.py instead.

This version incorrectly keeps chain A of 1NU8 (Diprotin A binds to chain B).
See AGENTS.md Trap #3 for the chain issue.

Phase 1 Structure Preparation Script (v1 — deprecated)
- Preprocess 1NU8 (DPP-IV) -> DPP4_clean.pdb
- Extract GHRH(1-29) from 7CZ5 or 7V9M -> GHRH_1-29.pdb
- Build YAD tripeptide -> YAD_tripeptide.pdb
"""

import warnings
warnings.warn("prepare_structures.py is deprecated. Use prepare_structures_v2.py instead. "
              "This version keeps chain A, but Diprotin A binds to chain B.", DeprecationWarning)

import os
import sys

# Use the installed PyMOL
pymol_path = "/home/scroll/miniforge3/lib/python3.12/site-packages"
if os.path.exists(pymol_path) and pymol_path not in sys.path:
    sys.path.insert(0, pymol_path)

try:
    import pymol
    from pymol import cmd
    pymol.finish_launching(["pymol", "-cq"])
except Exception as e:
    print(f"PyMOL import error: {e}")
    sys.exit(1)

BASE = "/home/scroll/personal/cjc-1295/workspace"
STRUCT = f"{BASE}/structures"
STEP0 = f"{BASE}/step0"
STEP1 = f"{BASE}/step1"

os.makedirs(STEP0, exist_ok=True)
os.makedirs(STEP1, exist_ok=True)

print("=" * 60)
print("Step 1.1: Preprocessing 1NU8 -> DPP4_clean.pdb")
print("=" * 60)

cmd.reinitialize()
cmd.load(f"{STRUCT}/1NU8.pdb", "1NU8")

# Remove solvent and unwanted heteroatoms
cmd.remove("solvent")
cmd.remove("resn HOH")
cmd.remove("resn SO4")
cmd.remove("resn NAG")

# Keep only chain A (one monomer) and chain D (Diprotin A reference)
# We'll save two versions: with and without Diprotin A
cmd.create("DPP4_with_ligand", "chain A or chain D")
cmd.create("DPP4_clean", "chain A")

# Save clean receptor
cmd.save(f"{STEP1}/DPP4_clean.pdb", "DPP4_clean")
# Save with ligand for reference
cmd.save(f"{STEP1}/DPP4_with_diprotinA.pdb", "DPP4_with_ligand")

# Report key residues
cmd.select("catalytic_triad", "resi 630+708+740 and chain A")
cmd.select("glu_motif", "resi 205+206 and chain A")
cmd.select("s1_pocket", "resi 631+656+659+662+666+711 and chain A")
cmd.select("oxy_hole", "resi 547+631 and chain A")

print(f"  Saved: {STEP1}/DPP4_clean.pdb")
print(f"  Saved: {STEP1}/DPP4_with_diprotinA.pdb")

# Check Diprotin A position relative to catalytic residues
print("\n--- Diprotin A (chain D) position check ---")
cmd.iterate_state(1, "chain D and name CA", "print('  DiprotinA res:',resn,resi)")

print("\n" + "=" * 60)
print("Step 1.2: Extract GHRH(1-29) from 7CZ5")
print("=" * 60)

cmd.reinitialize()
cmd.load(f"{STRUCT}/7CZ5.pdb", "7CZ5")

# Chain P is GHRH (residues 1-44 in the PDB correspond to GHRH 32-75 in UNP)
# We need GHRH(1-29) which corresponds to the first 29 residues of chain P
cmd.create("GHRH_full", "chain P")
cmd.create("GHRH_1-29", "chain P and resi 1-29")

# Save
cmd.save(f"{STEP1}/GHRH_1-29_from_7CZ5.pdb", "GHRH_1-29")

# Report sequence
print("  Extracted GHRH(1-29) from 7CZ5 chain P")
cmd.iterate_state(1, "GHRH_1-29 and name CA", "print('  ',resi,resn)")

print("\n" + "=" * 60)
print("Step 0: Build YAD tripeptide")
print("=" * 60)

cmd.reinitialize()
# CRITICAL FIX: Do NOT use cmd.fab("YAD", ss=0) which builds a linear
# peptide with phi=psi=0 (Ramachandran forbidden). Instead, extract
# the natural YAD conformation from the GHRH crystal structure (7CZ5).
# This preserves native backbone geometry (phi~96, psi~40 for Ala2)
# and prevents pair_fit() in prepare_docking_start.py from flattening
# the GHRH backbone.
cmd.load(f"{STEP1}/GHRH_1-29_from_7CZ5.pdb", "GHRH_ref")
cmd.create("YAD", "GHRH_ref and resi 1-3")
cmd.delete("GHRH_ref")

# Ensure protonation state: N-terminus NH3+, C-terminus COO-
# (the extracted YAD already has proper protonation from PDB)
cmd.save(f"{STEP0}/YAD_tripeptide.pdb", "YAD")
print(f"  Saved: {STEP0}/YAD_tripeptide.pdb")
print("  NOTE: YAD extracted from native GHRH(1-29) crystal structure")
print("        (avoids linear fab() artifact with phi=psi=0)")

# Also save a version aligned to Diprotin A in 1NU8 for quick validation
cmd.reinitialize()
cmd.load(f"{STEP1}/DPP4_with_diprotinA.pdb", "ref")
cmd.load(f"{STEP0}/YAD_tripeptide.pdb", "yad")

# Align YAD to Diprotin A (Ile-Pro-Ile)
# We align based on backbone atoms of the central residue (Pro2 of Diprotin A ~ Ala2 of YAD)
# This is a rough alignment for quick visualization
cmd.align("yad and name CA+C+N+O", "ref and chain D and name CA+C+N+O")

cmd.save(f"{STEP0}/YAD_aligned_to_DPP4.pdb", "yad")
cmd.save(f"{STEP0}/DPP4_with_YAD.pdb", "ref or yad")

print(f"  Saved: {STEP0}/YAD_aligned_to_DPP4.pdb")
print(f"  Saved: {STEP0}/DPP4_with_YAD.pdb")

print("\n" + "=" * 60)
print("Structure preparation complete!")
print("=" * 60)

cmd.quit()
