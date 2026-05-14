#!/usr/bin/env python3
"""
Prepare the starting structure for Rosetta FlexPepDock.
Align GHRH(1-29) N-terminus to the YAD tripeptide position in DPP-IV.

CRITICAL FIX (2026-05-14):
- Previous version used cmd.fab("YAD", ss=0) which builds a LINEAR peptide
  with phi=psi=0 (Ramachandran forbidden region).
- pair_fit() then forced GHRH backbone to match this linear reference,
  flattening the natural backbone and creating an unphysical starting
  conformation that caused GROMACS segfaults with dt=0.002.
- FIX: YAD_ref is now extracted from native GHRH(1-29) crystal structure
  (see prepare_structures.py). pair_fit() preserves natural backbone.
- Added phi/psi validation to catch any future backbone distortion.
"""

import os, sys, math

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

# -------------------------------------------------------------------------
# BACKBONE ALIGNMENT: pair_fit() rotates/translates GHRH so its N-terminal
# YAD matches the reference YAD position in DPP4 active site.
# IMPORTANT: This is a RIGID transformation. If YAD_ref has a natural
# backbone conformation (from crystal structure), GHRH backbone is preserved.
# If YAD_ref is linear (from fab(), ss=0), GHRH gets flattened.
# -------------------------------------------------------------------------
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

# -------------------------------------------------------------------------
# PHI/PSI VALIDATION: Check backbone dihedrals of GHRH residues 1-3
# to detect if pair_fit() distorted the backbone.
# Allowed regions (approximate for general amino acids):
#   Alpha-helix:  phi ~ -57, psi ~ -47
#   Beta-sheet:   phi ~ -120, psi ~ 120
#   Forbidden:    phi ~ 0, psi ~ 0  (linear artifact from fab())
# -------------------------------------------------------------------------
def get_atom_coords(obj, resi, name):
    """Return (x, y, z) for a specific atom."""
    model = cmd.get_model(f"{obj} and resi {resi} and name {name}")
    if model.atom:
        return model.atom[0].coord
    return None

def dihedral(p1, p2, p3, p4):
    """Calculate dihedral angle in degrees."""
    def vec_sub(a, b):
        return (a[0]-b[0], a[1]-b[1], a[2]-b[2])
    def vec_cross(a, b):
        return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
    def vec_dot(a, b):
        return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
    def vec_norm(a):
        return math.sqrt(a[0]**2 + a[1]**2 + a[2]**2)

    b1 = vec_sub(p2, p1)
    b2 = vec_sub(p3, p2)
    b3 = vec_sub(p4, p3)

    n1 = vec_cross(b1, b2)
    n2 = vec_cross(b2, b3)
    m1 = vec_cross(n1, b2)

    y = vec_norm(m1)
    x = vec_dot(n1, n2)
    angle = math.degrees(math.atan2(y, x))
    if vec_dot(n1, b3) > 0:
        angle = -angle
    return angle

print("\n--- GHRH N-terminus phi/psi validation ---")
phi_psi_ok = True
for resi in [2, 3]:
    prev_c = get_atom_coords("GHRH", resi-1, "C")
    n      = get_atom_coords("GHRH", resi,   "N")
    ca     = get_atom_coords("GHRH", resi,   "CA")
    c      = get_atom_coords("GHRH", resi,   "C")
    next_n = get_atom_coords("GHRH", resi+1, "N")

    if all(v is not None for v in [prev_c, n, ca, c, next_n]):
        phi = dihedral(prev_c, n, ca, c)
        psi = dihedral(n, ca, c, next_n)
        omega = dihedral(ca, c, next_n, get_atom_coords("GHRH", resi+1, "CA"))
        status = "OK"
        # Flag suspicious regions
        if abs(phi) < 30 and abs(psi) < 30:
            status = "WARNING: linear/fab artifact (phi~0, psi~0)"
            phi_psi_ok = False
        elif abs(omega) < 30 or abs(abs(omega) - 180) > 30:
            status = "WARNING: unusual omega"
            phi_psi_ok = False
        print(f"  Res {resi}: phi={phi:7.1f}°, psi={psi:7.1f}°, omega={omega:7.1f}°  [{status}]")
    else:
        print(f"  Res {resi}: Could not compute (missing atoms)")

if not phi_psi_ok:
    print("\n*** WARNING: Backbone dihedrals look suspicious.")
    print("*** If you see phi~0, psi~0, check that YAD_ref was NOT built with")
    print("*** cmd.fab('YAD', ss=0). Use native YAD from crystal structure instead.")
    print("*** See: prepare_structures.py for the correct YAD generation.")
else:
    print("  Backbone dihedrals look reasonable.")

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
