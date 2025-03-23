import re
import sys
from collections import defaultdict

"""
This script analyzes topology and forcefield files to identify unique angle and dihedral types in CHARMM36.
It extracts all angles and dihedrals from the topology and matches them to their corresponding types in the forcefield.
The script reports both successfully matched angle/dihedral types and any that could not be matched.

Warning: The logic in this script does not account for inconsistent naming. See `match_dihedral_to_dihedral_type`.
If needed, the fix for this would only require additional loops in the `x_to_x_type` functions.
"""

# =================================================================
FORCEFIELD_FILE = "charmm36.itp"  # Replace with your forcefield file path
TOPOLOGY_FILE = "popc.itp"      # Replace with your topology file path
# =================================================================

def parse_atom_types(topology_file):
    """Parse the atom types from the topology file"""
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
    """Parse the dihedraltypes from the forcefield file, including all parameters"""
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
                    i, j, k, l = parts[:4]
                    func = int(parts[4])
                    
                    # Collect all remaining parameters
                    params = parts[5:] if len(parts) > 5 else []
                    
                    # Store as a tuple: (atom_types, func, parameters)
                    dihedraltypes.append(((i, j, k, l), func, tuple(params)))
    
    return dihedraltypes

def parse_angletypes(forcefield_file):
    """Parse the angletypes from the forcefield file, including all parameters"""
    angletypes = []
    in_angletypes_section = False
    
    with open(forcefield_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            if line.startswith('[ angletypes ]'):
                in_angletypes_section = True
                continue
            
            if in_angletypes_section and line.startswith('['):
                in_angletypes_section = False
                continue
            
            if in_angletypes_section and line and not line.startswith(';'):
                parts = line.split()
                if len(parts) >= 4:
                    i, j, k = parts[:3]
                    func = int(parts[3])
                    
                    # Collect all remaining parameters
                    params = parts[4:] if len(parts) > 4 else []
                    
                    # Store as a tuple: (atom_types, func, parameters)
                    angletypes.append(((i, j, k), func, tuple(params)))
    
    return angletypes

def parse_dihedrals(topology_file):
    """Parse the dihedrals from the topology file"""
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

def parse_angles(topology_file):
    """Parse the angles from the topology file"""
    angles = []
    in_angles_section = False
    
    with open(topology_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            if line.startswith('[ angles ]'):
                in_angles_section = True
                continue
            
            if in_angles_section and line.startswith('['):
                in_angles_section = False
                continue
            
            if in_angles_section and line and not line.startswith(';'):
                parts = line.split()
                if len(parts) >= 4:
                    ai, aj, ak, func = map(int, parts[:4])
                    angles.append((ai, aj, ak, func))
    
    return angles

def match_dihedral_to_dihedral_type(dihedral_atoms, atom_types, dihedraltypes): 
    """
       Match a dihedral to its dihedral type. The logic in this function only accounts
       for topologies with consistent naming (i.e. if a dihedral is specified as ijkl, 
       jikl will be recognized as a unique torsion)
    """
    a1, a2, a3, a4 = [atom_types[a] for a in dihedral_atoms]
    
    # Find exact match
    for (i, j, k, l), func, params in dihedraltypes:
        if (i == a1 and j == a2 and k == a3 and l == a4) or \
           (i == a4 and j == a3 and k == a2 and l == a1):
            return ((i, j, k, l), func, params)
    
    # Find match with wildcards (X)
    for (i, j, k, l), func, params in dihedraltypes:
        if (j == a2 and k == a3) or (j == a3 and k == a2):
            if (i == a1 or i == "X") and (l == a4 or l == "X"):
                return ((i, j, k, l), func, params)
            if (i == a4 or i == "X") and (l == a1 or l == "X"):
                return ((i, j, k, l), func, params)
    
    # If no match found
    return None

def match_angle_to_angle_type(angle_atoms, atom_types, angletypes):
    """
       Match an angle to its angle type. The logic in this function only accounts
       for topologies with consistent naming (i.e. if an angle is specified as ijk, 
       jik will be recognized as a unique angle)
    """
    a1, a2, a3 = [atom_types[a] for a in angle_atoms]
    
    # Find exact match
    for (i, j, k), func, params in angletypes:
        if (i == a1 and j == a2 and k == a3) or \
           (i == a3 and j == a2 and k == a1):
            return ((i, j, k), func, params)
    
    # Find match with wildcards (X)
    for (i, j, k), func, params in angletypes:
        if j == a2:
            if (i == a1 or i == "X") and (k == a3 or k == "X"):
                return ((i, j, k), func, params)
            if (i == a3 or i == "X") and (k == a1 or k == "X"):
                return ((i, j, k), func, params)
    
    # If no match found
    return None

def analyze_dihedrals_and_angles(forcefield_file, topology_file):
    """Analyze the dihedrals and angles in a topology file"""
    try:
        atom_types = parse_atom_types(topology_file)
        dihedraltypes = parse_dihedraltypes(forcefield_file)
        angletypes = parse_angletypes(forcefield_file)
        dihedrals = parse_dihedrals(topology_file)
        angles = parse_angles(topology_file)
        
        # Analyze dihedrals - now considering both naming and parameters
        unique_dihedral_types = set()
        unmatched_dihedrals = []
        dihedral_param_groups = defaultdict(list)
        
        for ai, aj, ak, al, func in dihedrals:
            dihedral_match = match_dihedral_to_dihedral_type((ai, aj, ak, al), atom_types, dihedraltypes)
            if dihedral_match:
                atom_types_tuple, func, params = dihedral_match
                
                # Add to unique types set (by full specification)
                unique_dihedral_types.add((atom_types_tuple, func, params))
                
                # Group by parameter values
                param_key = (func, params)
                if atom_types_tuple not in dihedral_param_groups[param_key]:
                    dihedral_param_groups[param_key].append(atom_types_tuple)
            else:
                unmatched_dihedrals.append((ai, aj, ak, al, func))
        
        # Analyze angles - now considering both naming and parameters
        unique_angle_types = set()
        unmatched_angles = []
        angle_param_groups = defaultdict(list)
        
        for ai, aj, ak, func in angles:
            angle_match = match_angle_to_angle_type((ai, aj, ak), atom_types, angletypes)
            if angle_match:
                atom_types_tuple, func, params = angle_match
                
                # Add to unique types set (by full specification)
                unique_angle_types.add((atom_types_tuple, func, params))
                
                # Group by parameter values
                param_key = (func, params)
                if atom_types_tuple not in angle_param_groups[param_key]:
                    angle_param_groups[param_key].append(atom_types_tuple)
            else:
                unmatched_angles.append((ai, aj, ak, func))
        
        return unique_dihedral_types, unmatched_dihedrals, unique_angle_types, unmatched_angles, dihedral_param_groups, angle_param_groups
    except Exception as e:
        print(f"Error during analysis: {e}")
        return set(), [], set(), [], defaultdict(list), defaultdict(list)

def main():
    print(f"Using forcefield file: {FORCEFIELD_FILE}")
    print(f"Using topology file: {TOPOLOGY_FILE}")
    
    try:
        result = analyze_dihedrals_and_angles(FORCEFIELD_FILE, TOPOLOGY_FILE)
        unique_dihedral_types, unmatched_dihedrals, unique_angle_types, unmatched_angles, dihedral_param_groups, angle_param_groups = result
        
        # Print angle analysis results
        print("\n== ANGLE ANALYSIS ==")
        print(f"Found {len(unique_angle_types)} unique angle types (by atom naming):")
        for i, (atom_types, func, params) in enumerate(sorted(unique_angle_types), 1):
            a, b, c = atom_types
            param_str = " ".join(params)
            print(f"{i}. {a} {b} {c} (func {func}) {param_str}")
        
        print(f"\nFound {len(angle_param_groups)} unique angle parameter sets:")
        for i, ((func, params), atom_types_list) in enumerate(sorted(angle_param_groups.items()), 1):
            param_str = " ".join(params)
            print(f"{i}. Parameters: func {func} {param_str}")
            print(f"   Used by {len(atom_types_list)} different atom type combinations:")
            for j, (a, b, c) in enumerate(atom_types_list[:5], 1):
                print(f"   {j}. {a} {b} {c}")
            if len(atom_types_list) > 5:
                print(f"   ... and {len(atom_types_list) - 5} more")
        
        if unmatched_angles:
            print(f"\nWarning: {len(unmatched_angles)} angles could not be matched to an angle type.")
            print("First 5 unmatched angles:")
            for i, (ai, aj, ak, func) in enumerate(unmatched_angles[:5], 1):
                print(f"{i}. Atoms {ai}-{aj}-{ak} (func {func})")
        
        # Print dihedral analysis results
        print("\n== DIHEDRAL ANALYSIS ==")
        print(f"Found {len(unique_dihedral_types)} unique dihedral types (by atom naming):")
        for i, (atom_types, func, params) in enumerate(sorted(unique_dihedral_types), 1):
            a, b, c, d = atom_types
            param_str = " ".join(params)
            print(f"{i}. {a} {b} {c} {d} (func {func}) {param_str}")
        
        print(f"\nFound {len(dihedral_param_groups)} unique dihedral parameter sets:")
        for i, ((func, params), atom_types_list) in enumerate(sorted(dihedral_param_groups.items()), 1):
            param_str = " ".join(params)
            print(f"{i}. Parameters: func {func} {param_str}")
            print(f"   Used by {len(atom_types_list)} different atom type combinations:")
            for j, (a, b, c, d) in enumerate(atom_types_list[:5], 1):
                print(f"   {j}. {a} {b} {c} {d}")
            if len(atom_types_list) > 5:
                print(f"   ... and {len(atom_types_list) - 5} more")
        
        if unmatched_dihedrals:
            print(f"\nWarning: {len(unmatched_dihedrals)} dihedrals could not be matched to a dihedral type.")
            print("First 5 unmatched dihedrals:")
            for i, (ai, aj, ak, al, func) in enumerate(unmatched_dihedrals[:5], 1):
                print(f"{i}. Atoms {ai}-{aj}-{ak}-{al} (func {func})")
    
    # These exceptions may not print correctly ¯\_(ツ)_/¯
    except FileNotFoundError as e:
        print(f"Error: {e}")
        print("Please modify the file paths at the top of the script to point to valid files.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == "__main__":
    main()