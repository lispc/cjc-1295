#!/usr/bin/env python3
"""
Verify if OpenMM + openmmforcefields supports D-Ala (DAL).
Fix N-terminus by removing existing H and letting Modeller add correct ones.
"""

import sys
import os

sys.path.insert(0, "/home/scroll/miniforge3/lib/python3.13/site-packages")

try:
    from openmm import app
    import openmm
    from openmmforcefields.generators import SystemGenerator
    print(f"OpenMM version: {openmm.__version__}")
except ImportError as e:
    print(f"OpenMM not available yet: {e}")
    sys.exit(1)

def fix_and_test(pdb_path, label):
    print(f"\n{'='*60}")
    print(f"Testing: {label}")
    print(f"{'='*60}")
    
    pdb = app.PDBFile(pdb_path)
    print(f"Loaded {pdb.topology.getNumAtoms()} atoms")
    
    # Remove N-terminus H if present (let Modeller add correct NH3+)
    atoms_to_remove = []
    for res in pdb.topology.residues():
        if res.index == 0:  # N-terminus
            for atom in res.atoms():
                if atom.name == "H":
                    atoms_to_remove.append(atom)
                    print(f"  Removing existing N-term H from {res.name}")
    
    modeller = app.Modeller(pdb.topology, pdb.positions)
    if atoms_to_remove:
        modeller.delete(atoms_to_remove)
    
    # Test 1: amber14sb
    print("\n--- amber14sb ---")
    try:
        ff = app.ForceField("amber/protein.ff14SB.xml")
        modeller_amber = app.Modeller(modeller.topology, modeller.positions)
        modeller_amber.addHydrogens(forcefield=ff)
        system = ff.createSystem(modeller_amber.topology)
        print(f"  SUCCESS: {system.getNumParticles()} particles")
        return "amber14sb OK"
    except Exception as e:
        print(f"  FAILED: {e}")
    
    # Test 2: charmm36 with D-amino acid patch
    print("\n--- charmm36 + D-patch ---")
    try:
        ff_d = app.ForceField("charmm/charmm36_protein.xml", "charmm/charmm36_protein_d.xml")
        modeller_charmm = app.Modeller(modeller.topology, modeller.positions)
        modeller_charmm.addHydrogens(forcefield=ff_d)
        system = ff_d.createSystem(modeller_charmm.topology)
        print(f"  SUCCESS: {system.getNumParticles()} particles")
        return "charmm36 OK"
    except Exception as e:
        print(f"  FAILED: {e}")
    
    return "Both failed"

# Test WT
result_wt = fix_and_test(
    "/home/scroll/personal/cjc-1295/workspace/step1/GHRH_1-29.pdb",
    "WT GHRH(1-29)"
)

# Test D-Ala2
result_dala = fix_and_test(
    "/home/scroll/personal/cjc-1295/workspace/step1/GHRH_1-29_DAla2.pdb",
    "D-Ala2 Mutant"
)

print(f"\n{'='*60}")
print("Summary:")
print(f"  WT: {result_wt}")
print(f"  D-Ala2: {result_dala}")
print(f"{'='*60}")
