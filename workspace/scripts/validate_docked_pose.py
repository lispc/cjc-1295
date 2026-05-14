#!/usr/bin/env python3
"""
Validate the top-scoring Rosetta FlexPepDock pose against catalytic geometry criteria.
Compare with Step 0 YAD tripeptide benchmarks.

Chain-agnostic: looks up residues by residue number + residue name,
so it works regardless of Rosetta's chain assignment.
"""

import numpy as np
import sys

def read_pdb_atoms(path):
    """Read PDB and return dict keyed by (resi, resn, name). 
    Also returns chain_map: {(resi, resn): chain} for debugging."""
    atoms = {}
    chain_map = {}
    with open(path) as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                name = line[12:16].strip()
                resn = line[17:20].strip()
                chain = line[21]
                resi = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atoms[(resi, resn, name)] = np.array([x, y, z])
                chain_map[(resi, resn)] = chain
    return atoms, chain_map

def get_atom(atoms, resi, resn, name):
    """Fetch atom by (resi, resn, name). Supports wildcard resn for ALA/DALA."""
    key = (resi, resn, name)
    if key in atoms:
        return atoms[key]
    # Fallback: try DALA if ALA not found (or vice versa)
    alt_resn = "DALA" if resn == "ALA" else "ALA"
    alt_key = (resi, alt_resn, name)
    if alt_key in atoms:
        return atoms[alt_key]
    # Try any residue at this position with matching atom name
    for (ri, rn, an), pos in atoms.items():
        if ri == resi and an == name:
            return pos
    raise KeyError(f"Atom not found: residue {resi} {resn} atom {name}")

def dist(a, b):
    return np.linalg.norm(a - b)

def angle(a, b, c):
    ba = a - b
    bc = c - b
    cos = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    cos = np.clip(cos, -1.0, 1.0)
    return np.degrees(np.arccos(cos))

def main():
    pose_path = sys.argv[1] if len(sys.argv) > 1 else "/home/scroll/personal/cjc-1295/workspace/step2/prepacked_DPP4_GHRH_start_0001_0036_3.pdb"
    atoms, chain_map = read_pdb_atoms(pose_path)

    # DPP-IV catalytic residues (detect chains automatically for display)
    ser630_chain = chain_map.get((630, "SER"), "?")
    tyr547_chain = chain_map.get((547, "TYR"), "?")
    glu205_chain = chain_map.get((205, "GLU"), "?")
    glu206_chain = chain_map.get((206, "GLU"), "?")

    Ser630_OG = get_atom(atoms, 630, "SER", "OG")
    Tyr547_OH = get_atom(atoms, 547, "TYR", "OH")
    Glu205_OE1 = get_atom(atoms, 205, "GLU", "OE1")
    Glu205_OE2 = get_atom(atoms, 205, "GLU", "OE2")
    Glu206_OE1 = get_atom(atoms, 206, "GLU", "OE1")
    Glu206_OE2 = get_atom(atoms, 206, "GLU", "OE2")

    # GHRH N-terminus — support both ALA and DALA at position 2
    tyr1_chain = chain_map.get((1, "TYR"), "?")
    ala2_chain = chain_map.get((2, "ALA"), chain_map.get((2, "DALA"), "?"))

    Tyr1_N = get_atom(atoms, 1, "TYR", "N")
    Ala2_N = get_atom(atoms, 2, "ALA", "N")
    Ala2_CA = get_atom(atoms, 2, "ALA", "CA")
    Ala2_C = get_atom(atoms, 2, "ALA", "C")
    Ala2_O = get_atom(atoms, 2, "ALA", "O")

    # Measurements
    d_ser630_ala2c = dist(Ser630_OG, Ala2_C)
    d_ser630_ala2o = dist(Ser630_OG, Ala2_O)
    attack_angle = angle(Ser630_OG, Ala2_C, Ala2_N)

    glu_oe = [Glu205_OE1, Glu205_OE2, Glu206_OE1, Glu206_OE2]
    d_tyr1_glu = min(dist(Tyr1_N, oe) for oe in glu_oe)
    d_ala2o_tyr547 = dist(Ala2_O, Tyr547_OH)

    print(f"\n{'='*60}")
    print(f"  DOCKED POSE VALIDATION: {pose_path.split('/')[-1]}")
    print(f"{'='*60}")
    print(f"  Detected chains: Ser630=[{ser630_chain}], Tyr1=[{tyr1_chain}], Ala2=[{ala2_chain}]")

    print(f"\n  [Catalytic Geometry]")
    print(f"  Ser630 OG → Ala2 C=O C:  {d_ser630_ala2c:.2f} Å")
    print(f"  Ser630 OG → Ala2 C=O O:  {d_ser630_ala2o:.2f} Å")
    print(f"  Attack angle ∠(OG–C–N):   {attack_angle:.1f}°")

    print(f"\n  [Substrate Positioning]")
    print(f"  Tyr1 N → nearest Glu Oε: {d_tyr1_glu:.2f} Å")
    print(f"  Ala2 O → Tyr547 OH:      {d_ala2o_tyr547:.2f} Å")

    # Benchmarks from Step 0 YAD tripeptide
    print(f"\n  [Benchmark Comparison — Step 0 YAD Tripeptide]")
    print(f"  {'Metric':<35} {'Step0':>8} {'Docked':>8} {'Status':>10}")
    print(f"  {'-'*62}")

    def status(val, low, high, label):
        if low <= val <= high:
            return "PASS"
        return f"{label}"

    s1 = status(d_ser630_ala2c, 2.0, 4.0, "FAR")
    s2 = status(attack_angle, 80, 120, "BAD")
    s3 = status(d_tyr1_glu, 0, 3.5, "WEAK")
    s4 = status(d_ala2o_tyr547, 2.5, 4.0, "LOST")

    print(f"  {'Ser630 OG → Ala2 C (Å)':<35} {'2.53':>8} {d_ser630_ala2c:>8.2f} {s1:>10}")
    print(f"  {'Attack angle (°)':<35} {'113.6':>8} {attack_angle:>8.1f} {s2:>10}")
    print(f"  {'Tyr1 N → Glu Oε (Å)':<35} {'1.87':>8} {d_tyr1_glu:>8.2f} {s3:>10}")
    print(f"  {'Ala2 O → Tyr547 OH (Å)':<35} {'3.15':>8} {d_ala2o_tyr547:>8.2f} {s4:>10}")

    # Overall assessment
    passed = sum([
        2.0 <= d_ser630_ala2c <= 4.0,
        80 <= attack_angle <= 120,
        d_tyr1_glu <= 3.5,
        2.5 <= d_ala2o_tyr547 <= 4.0
    ])
    print(f"\n  {'='*62}")
    print(f"  PASS: {passed}/4 criteria")
    if passed >= 3:
        print(f"  Pose is CATALYTICALLY COMPETENT — ready for MD")
    else:
        print(f"  Pose may need refinement")
    print(f"{'='*60}\n")

    return passed

if __name__ == '__main__':
    main()
