import re
import sys
from collections import defaultdict

"""This script analyzes topology and forcefield files to identify unique angle and dihedral types for a POPC lipid in CHARMM36.
It extracts all angles and dihedrals from the topology and matches them to their corresponding types in the forcefield.
The script reports both successfully matched angle/dihedral types and any that could not be matched.

Usage:
1. Modify the `FF_file` and `topology` variables at the top of the script
2. run
"""

FF_file = "charmm36.itp"
topology = "popc.itp"

def parse_atom_types(topology_file):
    atom_types = {}
    in_atoms_section = False
    
    with open(topology_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            if line.startswith('[ atoms ]'):
                in_atoms_section = True
                continue
            
            if in_atoms_section and line.startswith('['):
                in_atoms_section = False
                continue
            
            if in_atoms_section and line and not line.startswith(';'):
                parts = line.split()
                if len(parts) >= 3:
                    atom_nr = int(parts[0])
                    atom_type = parts[1]
                    atom_types[atom_nr] = atom_type
    
    return atom_types

def parse_dihedraltypes(forcefield_file):
    dihedraltypes = []
    in_dihedraltypes_section = False
    
    with open(forcefield_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            if line.startswith('[ dihedraltypes ]'):
                in_dihedraltypes_section = True
                continue
            
            if in_dihedraltypes_section and line.startswith('['):
                in_dihedraltypes_section = False
                continue
            
            if in_dihedraltypes_section and line and not line.startswith(';'):
                parts = line.split()
                if len(parts) >= 5:
                    i, j, k, l, func = parts[:5]
                    dihedraltypes.append((i, j, k, l, func))
    
    return dihedraltypes

def parse_dihedrals(topology_file):
    dihedrals = []
    in_dihedrals_section = False
    
    with open(topology_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            if line.startswith('[ dihedrals ]'):
                in_dihedrals_section = True
                continue
            
            if in_dihedrals_section and line.startswith('['):
                in_dihedrals_section = False
                continue
            
            if in_dihedrals_section and line and not line.startswith(';'):
                parts = line.split()
                if len(parts) >= 5:
                    ai, aj, ak, al, func = map(int, parts[:5])
                    dihedrals.append((ai, aj, ak, al, func))
    
    return dihedrals

def match_dihedral_to_dihedral_type(dihedral_atoms, atom_types, dihedraltypes):
    a1, a2, a3, a4 = [atom_types[a] for a in dihedral_atoms]
    
    # Find exact match
    for i, j, k, l, func in dihedraltypes:
        if (i == a1 and j == a2 and k == a3 and l == a4) or \
           (i == a4 and j == a3 and k == a2 and l == a1):
            return (i, j, k, l, func)
    
    # Find match with wildcards (X)
    for i, j, k, l, func in dihedraltypes:
        if (j == a2 and k == a3) or (j == a3 and k == a2):
            if (i == a1 or i == "X") and (l == a4 or l == "X"):
                return (i, j, k, l, func)
            if (i == a4 or i == "X") and (l == a1 or l == "X"):
                return (i, j, k, l, func)
    
    # If no match found
    wrong = "WRONG"
    return wrong

def analyze_dihedrals(forcefield_file, topology_file):
    try:
        atom_types = parse_atom_types(topology_file)
        dihedraltypes = parse_dihedraltypes(forcefield_file)
        dihedrals = parse_dihedrals(topology_file)
        
        unique_dihedral_types = set()
        unmatched_dihedrals = []
        
        for ai, aj, ak, al, func in dihedrals:
            dihedral_match = match_dihedral_to_dihedral_type((ai, aj, ak, al), atom_types, dihedraltypes)
            if dihedral_match:
                unique_dihedral_types.add(dihedral_match)
            else:
                unmatched_dihedrals.append((ai, aj, ak, al, func))
        
        return unique_dihedral_types, unmatched_dihedrals
    except Exception as e:
        print(f"Error during analysis: {e}")
        return set(), []

def main(): 
    unique_dihedral_types, unmatched_dihedrals = analyze_dihedrals(FF_file, topology)
    
    print(f"\nFound {len(unique_dihedral_types)} unique dihedral types:")
    for i, (a, b, c, d, func) in enumerate(sorted(unique_dihedral_types), 1):
        print(f"{i}. {a} {b} {c} {d} (func {func})")
    
    if unmatched_dihedrals:
        print(f"\nWarning: {len(unmatched_dihedrals)} dihedrals could not be matched to a dihedral type.")
        print("First 5 unmatched dihedrals:")
        for i, (ai, aj, ak, al, func) in enumerate(unmatched_dihedrals[:5], 1):
            print(f"{i}. Atoms {ai}-{aj}-{ak}-{al} (func {func})")

if __name__ == "__main__":
    main()