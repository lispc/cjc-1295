#!/usr/bin/env python3
"""
Convert D-Ala2 PDB to Rosetta-compatible format.
Changes ALA B 2 -> DAL B 2, ATOM -> HETATM for that residue.
"""
import sys

infile = sys.argv[1]
outfile = sys.argv[2]

with open(infile) as f:
    lines = f.readlines()

out = []
for line in lines:
    if line.startswith(("ATOM", "HETATM")):
        resn = line[17:20].strip()
        chain = line[21]
        resi = line[22:26].strip()
        if resn == "ALA" and chain == "B" and resi == "2":
            # Change to DAL and HETATM
            line = "HETATM" + line[6:17] + "DAL" + line[20:]
    out.append(line)

with open(outfile, "w") as f:
    f.writelines(out)

print(f"Converted {infile} -> {outfile}")
print("ALA B 2 -> DAL B 2 (HETATM)")
