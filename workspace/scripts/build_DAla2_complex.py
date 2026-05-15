
# ⚠️ CHIRALITY WARNING (2026-05-15):
# This script assumes the INPUT is L-Ala and OUTPUT is D-Ala2.
# However, the original template GHRH_1-29_from_7CZ5.pdb already had D-Ala at position 2.
# When run on that template, the output is actually L-Ala.
# See docs/CHIRALITY_CORRECTION.md for details.

#!/usr/bin/env python3
"""
Build D-Ala2 mutant from WT DPP-IV/GHRH complex.
Mirrors CB/HA across the N-CA-C plane.
"""

import numpy as np
import sys

def read_pdb(path):
    lines = []
    atoms = {}
    with open(path) as f:
        for line in f:
            lines.append(line)
            if line.startswith(("ATOM", "HETATM")):
                serial = int(line[6:11])
                name = line[12:16].strip()
                resn = line[17:20].strip()
                chain = line[21]
                resi = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms[serial] = {
                    'name': name, 'resn': resn, 'resi': resi, 'chain': chain,
                    'x': x, 'y': y, 'z': z, 'line': line
                }
    return lines, atoms

def write_pdb(path, lines, modified_atoms):
    with open(path, 'w') as f:
        for line in lines:
            if line.startswith(("ATOM", "HETATM")):
                serial = int(line[6:11])
                if serial in modified_atoms:
                    a = modified_atoms[serial]
                    f.write(f"{line[:30]}{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}{line[54:]}")
                else:
                    f.write(line)
            else:
                f.write(line)

def point_to_plane_mirror(p, a, b, c):
    n = np.cross(b - a, c - a)
    n = n / np.linalg.norm(n)
    d = np.dot(p - a, n)
    return p - 2 * d * n

def dihedral(p1, p2, p3, p4):
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    m1 = n1 / np.linalg.norm(n1)
    m2 = n2 / np.linalg.norm(n2)
    y = np.dot(np.cross(m1, m2), b2 / np.linalg.norm(b2))
    x = np.dot(m1, m2)
    return np.degrees(np.arctan2(y, x))

def main():
    input_pdb = sys.argv[1] if len(sys.argv) > 1 else "/home/scroll/personal/cjc-1295/workspace/step3/em.pdb"
    output_pdb = sys.argv[2] if len(sys.argv) > 2 else "/home/scroll/personal/cjc-1295/workspace/step3/em_DAla2.pdb"
    
    lines, atoms = read_pdb(input_pdb)
    
    # Find Ala2 atoms (residue number 2, ALA)
    ala2_atoms = {a['name']: (a, serial) for serial, a in atoms.items() 
                  if a['resi'] == 2 and a['resn'] == 'ALA'}
    
    if len(ala2_atoms) < 9:
        print(f"Error: Found only {len(ala2_atoms)} Ala2 atoms, expected ~10")
        print("Available:", list(ala2_atoms.keys()))
        return 1
    
    N = np.array([ala2_atoms['N'][0]['x'], ala2_atoms['N'][0]['y'], ala2_atoms['N'][0]['z']])
    CA = np.array([ala2_atoms['CA'][0]['x'], ala2_atoms['CA'][0]['y'], ala2_atoms['CA'][0]['z']])
    C = np.array([ala2_atoms['C'][0]['x'], ala2_atoms['C'][0]['y'], ala2_atoms['C'][0]['z']])
    CB = np.array([ala2_atoms['CB'][0]['x'], ala2_atoms['CB'][0]['y'], ala2_atoms['CB'][0]['z']])
    HA = np.array([ala2_atoms['HA'][0]['x'], ala2_atoms['HA'][0]['y'], ala2_atoms['HA'][0]['z']])
    
    # Mirror CB and HA across N-CA-C plane
    CB_D = point_to_plane_mirror(CB, N, CA, C)
    HA_D = point_to_plane_mirror(HA, N, CA, C)
    
    # Verify dihedral flip
    wt_chi = dihedral(N, CA, CB, C)
    da_chi = dihedral(N, CA, CB_D, C)
    
    print(f"WT  N-CA-CB-C dihedral: {wt_chi:.2f}°")
    print(f"D-Ala N-CA-CB-C dihedral: {da_chi:.2f}°")
    print(f"Sign flipped: {'YES' if wt_chi * da_chi < 0 else 'NO'}")
    
    # Update atoms
    modified = {}
    for name, new_pos in [('CB', CB_D), ('HA', HA_D)]:
        a, serial = ala2_atoms[name]
        modified[serial] = {'x': new_pos[0], 'y': new_pos[1], 'z': new_pos[2]}
    
    # Also mirror the methyl hydrogens (HB1, HB2, HB3)
    for hb_name in ['HB1', 'HB2', 'HB3']:
        if hb_name in ala2_atoms:
            hb = np.array([ala2_atoms[hb_name][0]['x'], ala2_atoms[hb_name][0]['y'], ala2_atoms[hb_name][0]['z']])
            hb_D = point_to_plane_mirror(hb, N, CA, C)
            a, serial = ala2_atoms[hb_name]
            modified[serial] = {'x': hb_D[0], 'y': hb_D[1], 'z': hb_D[2]}
    
    write_pdb(output_pdb, lines, modified)
    print(f"\nD-Ala2 mutant saved to: {output_pdb}")
    
    # Estimate steric clash: distance from flipped CB to nearby pocket atoms
    print("\n--- Steric Clash Check (CB flipped position) ---")
    pocket_residues = [629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 641, 642, 643, 644, 645, 646, 647, 648, 649, 650, 651, 652, 653, 654, 655, 656, 657, 658, 659, 660, 661]
    clash_atoms = []
    for serial, a in atoms.items():
        if a['resi'] in pocket_residues and a['resn'] != 'ALA':
            pos = np.array([a['x'], a['y'], a['z']])
            d = np.linalg.norm(pos - CB_D)
            if d < 4.0:
                clash_atoms.append((a['name'], a['resn'], a['resi'], d))
    
    clash_atoms.sort(key=lambda x: x[3])
    print(f"Found {len(clash_atoms)} atoms within 4.0 Å of flipped CB:")
    for name, resn, resi, d in clash_atoms[:15]:
        marker = "⚠️ CLASH" if d < 2.5 else "close"
        print(f"  {resn}{resi} {name}: {d:.2f} Å {marker}")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
