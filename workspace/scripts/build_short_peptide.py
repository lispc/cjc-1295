#!/usr/bin/env python3
"""Build DPP-IV + GHRH(1-10) complex for short peptide MD validation."""
import sys
sys.path.insert(0, '/home/scroll/miniforge3/lib/python3.13/site-packages')
import pymol2

COMPLEX_PDB = "/home/scroll/personal/cjc-1295/workspace/step2/GHRH_DPP4_docked_best.pdb"
OUT_PDB = "/home/scroll/personal/cjc-1295/workspace/step3/short_peptide_complex.pdb"

with pymol2.PyMOL() as pymol:
    cmd = pymol.cmd
    cmd.load(COMPLEX_PDB, "complex")

    # Determine chain assignments from the complex PDB
    # DPP-IV is typically chain A, GHRH is chain B
    # Create separate selections
    cmd.select("dppiv", "chain A")
    cmd.select("ghrh_full", "chain B")
    cmd.select("ghrh_1_10", "chain B and resi 1-10")

    # Rename chains for clean output
    cmd.alter("dppiv", "chain='A'")
    cmd.alter("ghrh_1_10", "chain='B'")

    # Save as combined PDB
    cmd.save(OUT_PDB, "dppiv or ghrh_1_10")

    n_dppiv = cmd.count_atoms("dppiv")
    n_ghrh = cmd.count_atoms("ghrh_1_10")
    print(f"DPP-IV atoms: {n_dppiv}")
    print(f"GHRH 1-10 atoms: {n_ghrh}")
    print(f"Written to: {OUT_PDB}")
