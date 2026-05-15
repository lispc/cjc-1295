
# ⚠️ CHIRALITY NOTE (2026-05-15):
# This FEP topology was built with D-Ala2 at lambda=0 and L-Ala at lambda=1
# (opposite of intended L→D direction). Result needs sign flip.
# See docs/CHIRALITY_CORRECTION.md

#!/usr/bin/env python3
"""
Generate dual-topology hybrid ITP for L-Ala2 <-> D-Ala2 alchemical FEP.

Strategy:
  - Shared atoms (N, H, CA, HA, C, O): same in both states
  - L-specific atoms (CB, HB1, HB2, HB3): normal in state A, dummy in state B
  - D-specific atoms (CBD, HB1D, HB2D, HB3D): dummy in state A, normal in state B

This script reads a standard Amber/GROMACS .itp for chain B (GHRH),
finds residue 2 (ALA), and generates a hybrid version.

Output: hybrid_chain_B.itp + hybrid_coords.gro
"""

import sys
import copy


def parse_itp(path):
    """Parse a GROMACS .itp file into sections."""
    sections = {}
    current = None
    header_lines = []
    with open(path) as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith('[') and stripped.endswith(']'):
                current = stripped[1:-1].strip()
                if current not in sections:
                    sections[current] = []
                sections[current].append(line)
            elif current is None:
                header_lines.append(line)
            else:
                sections[current].append(line)
    return header_lines, sections


def find_residue_atoms(atoms_lines, target_resi):
    """Find atom indices belonging to a specific residue."""
    atom_nums = []
    atom_names = {}
    for line in atoms_lines:
        if line.strip().startswith(';') or not line.strip():
            continue
        parts = line.split()
        if len(parts) >= 5:
            try:
                resi = int(parts[2])
                if resi == target_resi:
                    atom_num = int(parts[0])
                    atom_name = parts[4]
                    atom_nums.append(atom_num)
                    atom_names[atom_num] = atom_name
            except ValueError:
                pass
    return atom_nums, atom_names


def make_dummy_type(atom_type):
    """Return dummy variant of atom type.
    
    Uses 'DU' - a generic dummy atom type with zero vdW and charge.
    Must be defined in the topology's [ atomtypes ] section.
    """
    return "DU"


def generate_hybrid_itp(l_itp_path, d_itp_path, out_itp_path, target_resi=2):
    """Generate hybrid topology from L and D state topologies."""
    
    l_header, l_sections = parse_itp(l_itp_path)
    d_header, d_sections = parse_itp(d_itp_path)
    
    # Find atoms to duplicate in target residue
    l_atoms = l_sections.get('atoms', [])
    l_res_atoms, l_res_names = find_residue_atoms(l_atoms, target_resi)
    
    d_atoms = d_sections.get('atoms', [])
    d_res_atoms, d_res_names = find_residue_atoms(d_atoms, target_resi)
    
    print(f"L-state residue {target_resi} atoms: {l_res_names}")
    print(f"D-state residue {target_resi} atoms: {d_res_names}")
    
    # Identify side-chain atoms (non-backbone)
    backbone_atoms = {'N', 'H', 'HN', 'CA', 'HA', 'C', 'O'}
    l_sc_atoms = [n for n, name in l_res_names.items() if name not in backbone_atoms]
    d_sc_atoms = [n for n, name in d_res_names.items() if name not in backbone_atoms]
    
    print(f"L side-chain atoms: {[l_res_names[n] for n in l_sc_atoms]}")
    print(f"D side-chain atoms: {[d_res_names[n] for n in d_sc_atoms]}")
    
    # Build atom index mapping for new hybrid topology
    # We need to add D-specific atoms as new entries
    max_atom = 0
    for line in l_atoms:
        if line.strip().startswith(';'):
            continue
        parts = line.split()
        if len(parts) >= 1:
            try:
                max_atom = max(max_atom, int(parts[0]))
            except ValueError:
                pass
    
    next_atom = max_atom + 1
    d_atom_map = {}  # old D atom number -> new hybrid atom number
    d_name_map = {}  # D atom name -> new hybrid atom number
    
    for old_num in d_sc_atoms:
        d_atom_map[old_num] = next_atom
        d_name_map[d_res_names[old_num]] = next_atom
        next_atom += 1
    
    print(f"New D-specific atoms will be numbered: {list(d_atom_map.values())}")
    
    # Build hybrid [ atoms ] section
    hybrid_atoms = []
    for line in l_atoms:
        if line.strip().startswith(';'):
            hybrid_atoms.append(line)
            continue
        parts = line.split()
        if len(parts) < 7:
            hybrid_atoms.append(line)
            continue
        try:
            atom_num = int(parts[0])
        except ValueError:
            hybrid_atoms.append(line)
            continue
        
        if atom_num in l_res_names:
            atom_name = l_res_names[atom_num]
            if atom_name in backbone_atoms:
                # Backbone atom - shared, no B state
                hybrid_atoms.append(line)
            else:
                # L-specific side chain atom: add B state as dummy
                # Format: nr type resnr residue atom cgnr charge mass typeB chargeB massB
                if len(parts) >= 11:
                    # Already has B state
                    hybrid_atoms.append(line)
                else:
                    # Add dummy B state
                    base = ' '.join(parts[:8])
                    mass = parts[7]
                    hybrid_atoms.append(f"{base} DU 0.0000 {mass}\n")
        else:
            hybrid_atoms.append(line)
    
    # Add D-specific atoms
    for line in d_atoms:
        if line.strip().startswith(';'):
            continue
        parts = line.split()
        if len(parts) < 7:
            continue
        try:
            atom_num = int(parts[0])
        except ValueError:
            continue
        
        if atom_num in d_atom_map:
            new_num = d_atom_map[atom_num]
            # D-specific atom: dummy in A state, normal in B state
            # nr type resnr residue atom cgnr charge mass typeB chargeB massB
            base_fields = parts[:8]
            new_line = f"{new_num:5d}  DU      {parts[2]:5s} {parts[3]:5s} {parts[4]:6s} {parts[5]:5s}  0.0000 {parts[7]:>10s}  {parts[1]:6s} {parts[6]:>8s} {parts[7]:>10s}\n"
            hybrid_atoms.append(new_line)
    
    # Build hybrid [ bonds ] section
    hybrid_bonds = []
    for line in l_sections.get('bonds', []):
        if line.strip().startswith(';'):
            hybrid_bonds.append(line)
            continue
        parts = line.split()
        if len(parts) < 3:
            hybrid_bonds.append(line)
            continue
        try:
            ai, aj = int(parts[0]), int(parts[1])
        except ValueError:
            hybrid_bonds.append(line)
            continue
        
        # Check if both atoms are in target residue
        ai_in = ai in l_res_names
        aj_in = aj in l_res_names
        
        if ai_in and aj_in:
            # Bond within target residue
            # Check if this bond involves side-chain atoms
            ai_name = l_res_names.get(ai, '')
            aj_name = l_res_names.get(aj, '')
            if ai_name in backbone_atoms and aj_name in backbone_atoms:
                hybrid_bonds.append(line)  # Backbone bond, shared
            elif ai_name in backbone_atoms:
                # Backbone to L-sidechain: add B state (dummy bond)
                base = ' '.join(parts[:5])
                hybrid_bonds.append(f"{base} 0.0 0.0\n")
            elif aj_name in backbone_atoms:
                base = ' '.join(parts[:5])
                hybrid_bonds.append(f"{base} 0.0 0.0\n")
            else:
                # Sidechain-sidechain: add B state (dummy bond)
                base = ' '.join(parts[:5])
                hybrid_bonds.append(f"{base} 0.0 0.0\n")
        else:
            hybrid_bonds.append(line)
    
    # Add D-specific bonds (from D topology)
    for line in d_sections.get('bonds', []):
        if line.strip().startswith(';'):
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            ai, aj = int(parts[0]), int(parts[1])
        except ValueError:
            continue
        
        ai_in = ai in d_res_names
        aj_in = aj in d_res_names
        
        if ai_in and aj_in:
            ai_name = d_res_names.get(ai, '')
            aj_name = d_res_names.get(aj, '')
            
            # Map D atom numbers to hybrid numbers
            new_ai = d_atom_map.get(ai, ai)
            new_aj = d_atom_map.get(aj, aj)
            
            # Determine if this is a backbone-to-DSC or SC-to-DSC bond
            if ai_name in backbone_atoms and aj_name in backbone_atoms:
                continue  # Already added from L topology
            
            # D-specific bond: dummy in A, normal in B
            funct = parts[2]
            b0 = parts[3] if len(parts) > 3 else ''
            kb = parts[4] if len(parts) > 4 else ''
            new_line = f"{new_ai:5d} {new_aj:5d}  {funct}  0.0 0.0  {b0} {kb}\n"
            hybrid_bonds.append(new_line)
    
    # For angles, dihedrals, pairs - similar logic but much more complex
    # For now, we take a shortcut: since L and D topologies are identical
    # except for CB/HB coordinates, the bond/angle/dihedral lists are the same.
    # We just need to handle the ones involving CB/HB.
    
    # Actually, let me check if L and D angle/dihedral lists are identical
    l_angles = l_sections.get('angles', [])
    d_angles = d_sections.get('angles', [])
    
    # Simple comparison: count non-comment lines
    l_angle_count = sum(1 for l in l_angles if not l.strip().startswith(';') and l.strip())
    d_angle_count = sum(1 for l in d_angles if not l.strip().startswith(';') and l.strip())
    print(f"L angles: {l_angle_count}, D angles: {d_angle_count}")
    
    # If identical, we can copy the angle list and just modify the ones involving CB/HB
    
    # Build hybrid [ angles ]
    hybrid_angles = []
    for line in l_angles:
        if line.strip().startswith(';'):
            hybrid_angles.append(line)
            continue
        parts = line.split()
        if len(parts) < 4:
            hybrid_angles.append(line)
            continue
        try:
            ai, aj, ak = int(parts[0]), int(parts[1]), int(parts[2])
        except ValueError:
            hybrid_angles.append(line)
            continue
        
        atoms_in = [a for a in (ai, aj, ak) if a in l_res_names]
        if atoms_in:
            names = [l_res_names.get(a, '') for a in (ai, aj, ak)]
            sc_count = sum(1 for n in names if n not in backbone_atoms)
            if sc_count > 0:
                # Angle involves side chain: add dummy B state
                base = ' '.join(parts[:6])
                hybrid_angles.append(f"{base} 0.0 0.0\n")
            else:
                hybrid_angles.append(line)  # Backbone only
        else:
            hybrid_angles.append(line)
    
    # Add D-specific angles
    for line in d_angles:
        if line.strip().startswith(';'):
            continue
        parts = line.split()
        if len(parts) < 4:
            continue
        try:
            ai, aj, ak = int(parts[0]), int(parts[1]), int(parts[2])
        except ValueError:
            continue
        
        atoms_in = [a for a in (ai, aj, ak) if a in d_res_names]
        if atoms_in:
            names = [d_res_names.get(a, '') for a in (ai, aj, ak)]
            sc_count = sum(1 for n in names if n not in backbone_atoms)
            if sc_count == 0:
                continue  # Backbone angle, already added
            
            new_ai = d_atom_map.get(ai, ai)
            new_aj = d_atom_map.get(aj, aj)
            new_ak = d_atom_map.get(ak, ak)
            
            funct = parts[3]
            th0 = parts[4] if len(parts) > 4 else ''
            cth = parts[5] if len(parts) > 5 else ''
            new_line = f"{new_ai:5d} {new_aj:5d} {new_ak:5d}  {funct}  0.0 0.0  {th0} {cth}\n"
            hybrid_angles.append(new_line)
    
    # [ dihedrals ] - similar approach
    hybrid_dihedrals = []
    for line in l_sections.get('dihedrals', []):
        if line.strip().startswith(';'):
            hybrid_dihedrals.append(line)
            continue
        parts = line.split()
        if len(parts) < 5:
            hybrid_dihedrals.append(line)
            continue
        try:
            ai, aj, ak, al = int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3])
        except ValueError:
            hybrid_dihedrals.append(line)
            continue
        
        atoms_in = [a for a in (ai, aj, ak, al) if a in l_res_names]
        if atoms_in:
            names = [l_res_names.get(a, '') for a in (ai, aj, ak, al)]
            sc_count = sum(1 for n in names if n not in backbone_atoms)
            if sc_count > 0:
                base = ' '.join(parts[:8])  # funct + 4 params
                hybrid_dihedrals.append(f"{base} 0.0 0.0 0.0 0.0\n")
            else:
                hybrid_dihedrals.append(line)
        else:
            hybrid_dihedrals.append(line)
    
    # Add D-specific dihedrals
    for line in d_sections.get('dihedrals', []):
        if line.strip().startswith(';'):
            continue
        parts = line.split()
        if len(parts) < 5:
            continue
        try:
            ai, aj, ak, al = int(parts[0]), int(parts[1]), int(parts[2]), int(parts[3])
        except ValueError:
            continue
        
        atoms_in = [a for a in (ai, aj, ak, al) if a in d_res_names]
        if atoms_in:
            names = [d_res_names.get(a, '') for a in (ai, aj, ak, al)]
            sc_count = sum(1 for n in names if n not in backbone_atoms)
            if sc_count == 0:
                continue
            
            new_ai = d_atom_map.get(ai, ai)
            new_aj = d_atom_map.get(aj, aj)
            new_ak = d_atom_map.get(ak, ak)
            new_al = d_atom_map.get(al, al)
            
            funct = parts[4]
            params = ' '.join(parts[5:9]) if len(parts) > 5 else ''
            params_b = ' '.join(parts[5:9]) if len(parts) > 5 else ''
            new_line = f"{new_ai:5d} {new_aj:5d} {new_ak:5d} {new_al:5d}  {funct}  0.0 0.0 0.0 0.0  {params}\n"
            hybrid_dihedrals.append(new_line)
    
    # [ pairs ] - 1-4 interactions
    hybrid_pairs = []
    for line in l_sections.get('pairs', []):
        if line.strip().startswith(';'):
            hybrid_pairs.append(line)
            continue
        parts = line.split()
        if len(parts) < 3:
            hybrid_pairs.append(line)
            continue
        try:
            ai, aj = int(parts[0]), int(parts[1])
        except ValueError:
            hybrid_pairs.append(line)
            continue
        
        ai_in = ai in l_res_names
        aj_in = aj in l_res_names
        if ai_in and aj_in:
            names = [l_res_names.get(ai, ''), l_res_names.get(aj, '')]
            sc_count = sum(1 for n in names if n not in backbone_atoms)
            if sc_count > 0:
                base = ' '.join(parts[:4])
                hybrid_pairs.append(f"{base} 0.0 0.0\n")
            else:
                hybrid_pairs.append(line)
        else:
            hybrid_pairs.append(line)
    
    # Add D-specific pairs
    for line in d_sections.get('pairs', []):
        if line.strip().startswith(';'):
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        try:
            ai, aj = int(parts[0]), int(parts[1])
        except ValueError:
            continue
        
        ai_in = ai in d_res_names
        aj_in = aj in d_res_names
        if ai_in and aj_in:
            names = [d_res_names.get(ai, ''), d_res_names.get(aj, '')]
            sc_count = sum(1 for n in names if n not in backbone_atoms)
            if sc_count == 0:
                continue
            
            new_ai = d_atom_map.get(ai, ai)
            new_aj = d_atom_map.get(aj, aj)
            
            funct = parts[2]
            params = ' '.join(parts[3:5]) if len(parts) > 3 else ''
            new_line = f"{new_ai:5d} {new_aj:5d}  {funct}  0.0 0.0  {params}\n"
            hybrid_pairs.append(new_line)
    
    # [ exclusions ] - need to add exclusions between L-specific and D-specific atoms
    # to prevent them from interacting
    hybrid_exclusions = list(l_sections.get('exclusions', []))
    # Add exclusion entries for each L-SC atom with each D-SC atom
    for l_atom in l_sc_atoms:
        for d_atom in d_sc_atoms:
            new_d = d_atom_map[d_atom]
            hybrid_exclusions.append(f"{l_atom:5d} {new_d:5d}\n")
    
    # Write output
    with open(out_itp_path, 'w') as f:
        f.write("; Hybrid topology for L-Ala2 <-> D-Ala2 alchemical FEP\n")
        f.write("; Generated by generate_fep_topology.py\n")
        f.write(";\n")
        f.writelines(l_header)
        
        for section_name in l_sections:
            if section_name == 'atoms':
                f.writelines(hybrid_atoms)
            elif section_name == 'bonds':
                f.writelines(hybrid_bonds)
            elif section_name == 'angles':
                f.writelines(hybrid_angles)
            elif section_name == 'dihedrals':
                f.writelines(hybrid_dihedrals)
            elif section_name == 'pairs':
                f.writelines(hybrid_pairs)
            elif section_name == 'exclusions':
                f.writelines(hybrid_exclusions)
            else:
                f.writelines(l_sections[section_name])
    
    print(f"\nWrote hybrid topology: {out_itp_path}")
    print(f"Total atoms: {next_atom - 1}")
    print(f"Added {len(d_sc_atoms)} D-specific atoms")
    
    return d_atom_map


def generate_hybrid_gro(l_gro_path, d_gro_path, out_gro_path, d_atom_map, target_resi=2):
    """Generate hybrid GRO by appending D-specific atom coordinates."""
    
    with open(l_gro_path) as f:
        l_lines = f.readlines()
    
    with open(d_gro_path) as f:
        d_lines = f.readlines()
    
    # Parse L gro
    l_title = l_lines[0]
    l_natoms = int(l_lines[1].strip())
    l_coords = l_lines[2:2+l_natoms]
    l_box = l_lines[2+l_natoms]
    
    # Parse D gro
    d_natoms = int(d_lines[1].strip())
    d_coords = d_lines[2:2+d_natoms]
    
    # Find D-specific atoms in residue 2
    d_sc_coords = {}
    backbone_atoms = {'N', 'H', 'HN', 'CA', 'HA', 'C', 'O'}
    for line in d_coords:
        if len(line) < 20:
            continue
        resi_str = line[0:5].strip()
        try:
            resi = int(resi_str)
        except ValueError:
            continue
        if resi == target_resi:
            atom_name = line[10:15].strip()
            if atom_name not in backbone_atoms:
                # This is a D-specific side chain atom
                d_sc_coords[atom_name] = line
    
    # Build hybrid gro
    new_natoms = l_natoms + len(d_sc_coords)
    
    with open(out_gro_path, 'w') as f:
        f.write(l_title)
        f.write(f"{new_natoms:5d}\n")
        
        # Write all L atoms
        for line in l_coords:
            f.write(line)
        
        # Write D-specific atoms with new numbers
        atom_counter = l_natoms + 1
        for atom_name, line in sorted(d_sc_coords.items()):
            # Replace atom number and residue name (keep ALA)
            new_line = f"{atom_counter:5d}{line[5:10]}D{line[10:15]}{atom_counter:5d}{line[20:]}"
            f.write(new_line)
            atom_counter += 1
        
        f.write(l_box)
    
    print(f"Wrote hybrid GRO: {out_gro_path}")
    print(f"Added {len(d_sc_coords)} D-specific atoms to coordinates")


if __name__ == "__main__":
    if len(sys.argv) < 5:
        print("Usage: python generate_fep_topology.py <L_chain_B.itp> <D_chain_B.itp> <L.gro> <D.gro> [out_prefix]")
        sys.exit(1)
    
    l_itp = sys.argv[1]
    d_itp = sys.argv[2]
    l_gro = sys.argv[3]
    d_gro = sys.argv[4]
    prefix = sys.argv[5] if len(sys.argv) > 5 else "hybrid"
    
    d_atom_map = generate_hybrid_itp(l_itp, d_itp, f"{prefix}_chain_B.itp")
    generate_hybrid_gro(l_gro, d_gro, f"{prefix}.gro", d_atom_map)