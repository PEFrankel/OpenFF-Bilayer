

def read_atoms_section(filename):
    atoms = {}
    
    with open(filename, 'r') as file:
        lines = file.readlines()
        in_atoms_section = False
        
        for line in lines:
            if line.strip() == '[ atoms ]':
                in_atoms_section = True
                continue
            
            if in_atoms_section:
                if line.startswith('['): #skip comments or stop at header
                    break
                if line.startswith(';'):
                    continue
                
                parts = line.split()
                
                if len(parts) >= 5: # top parsing
                    index = int(parts[0])
                    atom_name = parts[4]
                    
                    if atom_name.startswith('C') or atom_name.startswith('H'): # store index if carbon or hydrogen
                        atoms[atom_name] = {'index': index, 'bonds': 0, 'mass': float(parts[7])}

    return atoms

def read_bonds_section(filename, atoms):
    with open(filename, 'r') as file:
        lines = file.readlines()
        in_bonds_section = False
        
        for line in lines: # find [ bonds ]
            if line.strip() == '[ bonds ]':
                in_bonds_section = True
                continue

            #in [ bonds ]
            if in_bonds_section:
                if line.startswith('['):
                    break
                if line.startswith(';'):
                    continue
                
                parts = line.split()
                if len(parts) >= 2:
                    index_a = int(parts[0])
                    index_b = int(parts[1])
                    
                    # Count h-bonds for each carbon
                    for atom_name, atom_info in atoms.items():
                        if atom_info['index'] == index_a and atom_name.startswith('C'):
                            if any(atom_info['index'] == index_b for atom_info in atoms.values() if atom_info['index'] == index_b and atom_info['index'] in [h_info['index'] for h_name, h_info in atoms.items() if h_name.startswith('H')]):
                                atoms[atom_name]['bonds'] += 1
                        elif atom_info['index'] == index_b and atom_name.startswith('C'):
                            if any(atom_info['index'] == index_a for atom_info in atoms.values() if atom_info['index'] == index_a and atom_info['index'] in [h_info['index'] for h_name, h_info in atoms.items() if h_name.startswith('H')]):
                                atoms[atom_name]['bonds'] += 1

def repartition_hydrogen_masses(atoms, original_H_mass, new_H_mass):
    H_difference = new_H_mass - original_H_mass

    for atom_name, atom_info in atoms.items():
        if atom_name.startswith('C'):
            bonded_hydrogens = atom_info['bonds']
            carbon_mass = atom_info['mass']
            new_carbon_mass = carbon_mass - (H_difference * bonded_hydrogens)
            atom_info['mass'] = new_carbon_mass
        
        elif atom_name.startswith('H'):
            atom_info['mass'] = new_H_mass

def main(filename):
    original_H_mass = 1.0080
    new_H_mass = 3.000000000000

    repartition_hydrogen_masses(atoms, original_H_mass, new_H_mass)

    atoms = read_atoms_section(filename)
    read_bonds_section(filename, atoms)

    for atom_name, atom_info in atoms.items():
        if atom_name.startswith('C'):
            print(f"{atom_name}: {atom_info['bonds']} bonds, new mass: {atom_info['mass']:.12f}")
        elif atom_name.startswith('H'):
            print(f"{atom_name}: new mass: {atom_info['mass']:.12f}")

filename = 'toppar/POPC.itp'
main(filename)
