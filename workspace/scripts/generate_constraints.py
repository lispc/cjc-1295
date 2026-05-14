#!/usr/bin/env python3
"""
Generate Rosetta CoordinateConstraint file for GHRH N-terminal CA atoms.
Usage: python generate_constraints.py <reference.pdb> <output.cst> [resi_list] [sd]
"""
import sys

def read_ca_coords(pdb_path, chain, resi_list):
    coords = {}
    with open(pdb_path) as f:
        for line in f:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            atom_name = line[12:16].strip()
            resn = line[17:20].strip()
            ch = line[21]
            resi = line[22:26].strip()
            if atom_name == "CA" and ch == chain and resi in resi_list:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords[resi] = (resn, x, y, z)
    return coords


if __name__ == "__main__":
    pdb_path = sys.argv[1]
    out_path = sys.argv[2]
    resi_list = sys.argv[3].split(",") if len(sys.argv) > 3 else ["1", "2", "3", "4", "5"]
    sd = sys.argv[4] if len(sys.argv) > 4 else "0.5"

    coords = read_ca_coords(pdb_path, "B", resi_list)

    with open(out_path, "w") as f:
        f.write("# CoordinateConstraints for GHRH N-terminal CA atoms\n")
        f.write(f"# Reference: {pdb_path}\n")
        f.write(f"# Chain B, residues {','.join(resi_list)}, sd={sd} Å\n\n")
        for resi in resi_list:
            if resi not in coords:
                print(f"WARNING: CA B {resi} not found in {pdb_path}")
                continue
            resn, x, y, z = coords[resi]
            # Rosetta CoordinateConstraint format:
            # CoordinateConstraint atom1_name res1_num atom2_name res2_num x y z func_type param
            # Using same atom as reference (self-reference with absolute coords)
            f.write(f"CoordinateConstraint CA {resi}B CA {resi}B {x:.3f} {y:.3f} {z:.3f} HARMONIC 0.0 {sd}\n")

    print(f"Wrote {len(coords)} constraints to {out_path}")
