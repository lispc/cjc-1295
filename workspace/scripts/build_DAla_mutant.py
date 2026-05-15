
# ⚠️ CHIRALITY WARNING (2026-05-15):
# This script assumes the INPUT is L-Ala and OUTPUT is D-Ala2.
# However, the original template GHRH_1-29_from_7CZ5.pdb already had D-Ala at position 2.
# When run on that template, the output is actually L-Ala.
# See docs/CHIRALITY_CORRECTION.md for details.

#!/usr/bin/env python3
"""
Build GHRH(1-29) D-Ala2 mutant by inverting chirality at residue 2.
Method: Mirror the CB atom across the N-CA-C plane.
"""

import os, sys, math

pymol_path = "/home/scroll/miniforge3/lib/python3.12/site-packages"
if os.path.exists(pymol_path) and pymol_path not in sys.path:
    sys.path.insert(0, pymol_path)

import pymol
from pymol import cmd
pymol.finish_launching(["pymol", "-cq"])

BASE = "/home/scroll/personal/cjc-1295/workspace"
STEP1 = f"{BASE}/step1"

cmd.reinitialize()
cmd.load(f"{STEP1}/GHRH_1-29.pdb", "GHRH")

# Manually mirror CB across the N-CA-C plane to invert chirality at CA
cmd.reinitialize()
cmd.load(f"{STEP1}/GHRH_1-29.pdb", "GHRH")

# Get coordinates
def get_atom_coords(obj, resi, name):
    model = cmd.get_model(f"{obj} and resi {resi} and name {name}")
    if model.atom:
        return model.atom[0].coord
    return None

N = get_atom_coords("GHRH", 2, "N")
CA = get_atom_coords("GHRH", 2, "CA")
C = get_atom_coords("GHRH", 2, "C")
CB = get_atom_coords("GHRH", 2, "CB")

print(f"Original N:  {N}")
print(f"Original CA: {CA}")
print(f"Original C:  {C}")
print(f"Original CB: {CB}")

# Mirror CB across the N-CA-C plane
# Plane normal = cross(N-CA, C-CA)
def vec_sub(a, b):
    return [a[i]-b[i] for i in range(3)]

def vec_cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]

def vec_dot(a, b):
    return sum(a[i]*b[i] for i in range(3))

def vec_norm(a):
    return math.sqrt(sum(x*x for x in a))

def vec_normalize(a):
    n = vec_norm(a)
    return [x/n for x in a]

v1 = vec_sub(N, CA)
v2 = vec_sub(C, CA)
normal = vec_cross(v1, v2)
normal = vec_normalize(normal)

# Distance from CB to plane
cb_vec = vec_sub(CB, CA)
dist = vec_dot(cb_vec, normal)

# Mirrored CB = CB - 2 * dist * normal
CB_new = [CB[i] - 2 * dist * normal[i] for i in range(3)]

print(f"Mirrored CB: {CB_new}")

# Update CB coordinates in PyMOL
cmd.alter_state(1, "GHRH and resi 2 and name CB", f"x={CB_new[0]}; y={CB_new[1]}; z={CB_new[2]}")

# Save the mutant
cmd.save(f"{STEP1}/GHRH_1-29_DAla2.pdb", "GHRH")
print(f"Saved: {STEP1}/GHRH_1-29_DAla2.pdb")

# Verify by computing N-CA-CB-C dihedral
def dihedral(p1, p2, p3, p4):
    b1 = vec_sub(p2, p1)
    b2 = vec_sub(p3, p2)
    b3 = vec_sub(p4, p3)
    
    b2n = vec_normalize(b2)
    
    v = vec_sub(b1, [vec_dot(b1, b2n) * b2n[i] for i in range(3)])
    w = vec_sub(b3, [vec_dot(b3, b2n) * b2n[i] for i in range(3)])
    
    v = vec_normalize(v)
    w = vec_normalize(w)
    
    x = vec_dot(v, w)
    y = vec_dot(vec_cross(b2n, v), w)
    
    return math.degrees(math.atan2(y, x))

# For L-Ala, N-CA-CB-C dihedral is typically around +60° (gauche+)
# For D-Ala, it should be around -60° (gauche-)
print(f"\nDihedral N-CA-CB-C (original): {dihedral(N, CA, CB, C):.1f}°")
print(f"Dihedral N-CA-CB-C (mutant):   {dihedral(N, CA, CB_new, C):.1f}°")

cmd.quit()
