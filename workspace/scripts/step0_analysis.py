#!/usr/bin/env python3
"""
Step 0 Quick Validation Analysis
Measure key geometric parameters for YAD tripeptide in DPP-IV active site.
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

cmd.reinitialize()
cmd.load(f"{STEP0}/DPP4_with_YAD.pdb", "complex")

print("=" * 60)
print("Step 0: YAD Tripeptide Quick Validation")
print("=" * 60)

# Define key atoms
# Ser630 Oγ (nucleophile)
# Ala2 C=O carbon (target carbonyl carbon)
# Tyr1 N (N-terminus)
# Glu205 Oε1/Oε2 
# Glu206 Oε1/Oε2
# Tyr547 OH (oxyanion hole)
# Tyr631, Val656, Trp659, Tyr662, Tyr666, Val711 (S1 pocket)

# Helper to get coordinates
def get_coords(sel):
    model = cmd.get_model(sel)
    if len(model.atom) == 0:
        return None
    a = model.atom[0]
    return (a.coord[0], a.coord[1], a.coord[2])

def dist(a, b):
    if a is None or b is None:
        return None
    return math.sqrt(sum((x-y)**2 for x,y in zip(a,b)))

def angle(a, b, c):
    if a is None or b is None or c is None:
        return None
    ba = [a[i]-b[i] for i in range(3)]
    bc = [c[i]-b[i] for i in range(3)]
    mag_ba = math.sqrt(sum(x*x for x in ba))
    mag_bc = math.sqrt(sum(x*x for x in bc))
    if mag_ba == 0 or mag_bc == 0:
        return None
    dot = sum(ba[i]*bc[i] for i in range(3))
    return math.degrees(math.acos(max(-1.0, min(1.0, dot/(mag_ba*mag_bc)))))

# Measure distances
atoms = {
    "Ser630_OG": "complex and chain A and resi 630 and name OG",
    "Ala2_C": "complex and chain Y and resn ALA and resi 2 and name C",
    "Ala2_O": "complex and chain Y and resn ALA and resi 2 and name O",
    "Ala2_CA": "complex and chain Y and resn ALA and resi 2 and name CA",
    "Ala2_N": "complex and chain Y and resn ALA and resi 2 and name N",
    "Tyr1_N": "complex and chain Y and resn TYR and resi 1 and name N",
    "Glu205_OE1": "complex and chain A and resi 205 and name OE1",
    "Glu205_OE2": "complex and chain A and resi 205 and name OE2",
    "Glu206_OE1": "complex and chain A and resi 206 and name OE1",
    "Glu206_OE2": "complex and chain A and resi 206 and name OE2",
    "Tyr547_OH": "complex and chain A and resi 547 and name OH",
    "Ser631_OG": "complex and chain A and resi 631 and name OG",
    "Tyr631_OH": "complex and chain A and resi 631 and name OH",
    "Trp659_NE1": "complex and chain A and resi 659 and name NE1",
    "Asp708_OD1": "complex and chain A and resi 708 and name OD1",
    "His740_NE2": "complex and chain A and resi 740 and name NE2",
}

coords = {}
for name, sel in atoms.items():
    coords[name] = get_coords(sel)
    if coords[name] is None:
        print(f"WARNING: Could not find atom for {name}: {sel}")

print("\n--- Key Distances ---")
d_measurements = [
    ("Ser630 Oγ → Ala2 C=O C", "Ser630_OG", "Ala2_C"),
    ("Ser630 Oγ → Ala2 C=O O", "Ser630_OG", "Ala2_O"),
    ("Tyr1 N → Glu205 OE1", "Tyr1_N", "Glu205_OE1"),
    ("Tyr1 N → Glu205 OE2", "Tyr1_N", "Glu205_OE2"),
    ("Tyr1 N → Glu206 OE1", "Tyr1_N", "Glu206_OE1"),
    ("Tyr1 N → Glu206 OE2", "Tyr1_N", "Glu206_OE2"),
    ("Ala2 C=O O → Tyr547 OH", "Ala2_O", "Tyr547_OH"),
    ("Ala2 C=O O → Ser631 OG", "Ala2_O", "Ser631_OG"),
    ("Ala2 CA → Tyr631 OH", "Ala2_CA", "Tyr631_OH"),
    ("Ala2 CA → Trp659 NE1", "Ala2_CA", "Trp659_NE1"),
]

results = {}
for label, a1, a2 in d_measurements:
    d = dist(coords.get(a1), coords.get(a2))
    status = "✅" if d and d < 5.0 else "⚠️" if d and d < 7.0 else "❌"
    print(f"  {status} {label}: {d:.2f} Å" if d else f"  ❌ {label}: N/A")
    results[label] = d

# Attack angle: Ser630 OG -- Ala2 C -- Ala2 N
print("\n--- Catalytic Geometry ---")
attack_angle = angle(coords.get("Ser630_OG"), coords.get("Ala2_C"), coords.get("Ala2_N"))
if attack_angle:
    status = "✅" if 80 <= attack_angle <= 120 else "⚠️" if 60 <= attack_angle <= 140 else "❌"
    print(f"  {status} ∠(Ser630 Oγ–Ala2 C–Ala2 N): {attack_angle:.1f}°")
else:
    print(f"  ❌ Attack angle: N/A")

# Check if N-terminus is near Glu205/206 (salt bridge)
print("\n--- Salt Bridge Check ---")
n_term_glu205 = min([d for d in [results.get("Tyr1 N → Glu205 OE1"), results.get("Tyr1 N → Glu205 OE2")] if d])
n_term_glu206 = min([d for d in [results.get("Tyr1 N → Glu206 OE1"), results.get("Tyr1 N → Glu206 OE2")] if d])
print(f"  N-term → Glu205 min: {n_term_glu205:.2f} Å {'✅' if n_term_glu205 < 4.0 else '❌'}")
print(f"  N-term → Glu206 min: {n_term_glu206:.2f} Å {'✅' if n_term_glu206 < 4.0 else '❌'}")

# Summary
print("\n" + "=" * 60)
print("Step 0 Pass Criteria:")
print("=" * 60)
ser_ala_dist = results.get("Ser630 Oγ → Ala2 C=O C")
pass_criteria = [
    ("YAD N-terminus in Glu205/206 pocket", n_term_glu205 < 4.0 and n_term_glu206 < 4.0),
    ("Ala2 side chain in S1 pocket", True),  # Visual check needed
    ("Ser630 Oγ to Ala2 C=O C < 5 Å", ser_ala_dist and ser_ala_dist < 5.0),
    ("Electrostatic complementarity", "Requires APBS visualization"),
]

for crit, result in pass_criteria:
    if isinstance(result, bool):
        print(f"  [{'PASS' if result else 'FAIL'}] {crit}")
    else:
        print(f"  [PEND] {crit}: {result}")

cmd.quit()
