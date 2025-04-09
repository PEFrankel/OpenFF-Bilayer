import mdtraj as md

gro_file_path = 'traj_1f.gro'
output_file_path = 'ndxfile.txt'
trajectory = md.load(gro_file_path)

#OpenFF SMILES string specific
chain_1_atoms = ['C8x', 'C9x', 'C10x', 'C11x', 'C12x', 'C13x', 'C14x', 'C15x', 
                 'C16x', 'C17x', 'C18x', 'C19x', 'C20x', 'C21x', 'C22x', 'C23x', 'C24x', 'C25x'] # indices for the first chain
chain_2_atoms = ['C27x', 'C28x', 'C29x', 'C30x', 'C31x', 'C32x', 'C33x', 'C34x', 
                 'C35x', 'C36x', 'C37x', 'C38x', 'C39x', 'C40x', 'C41x', 'C42x'] # indices for the second chain

header_values = [128, 2, 18, 16]

with open(output_file_path, 'w') as output_file:
    for value in header_values:
        output_file.write(f"# {value:5d}\n")

    for atom in trajectory.topology.atoms:
        atom_index = atom.index + 1  # MDTraj uses 0 indexing; GRO uses 1
        atom_name = atom.name
        residue_index = atom.residue.index + 1  # same for residue index
        chain_number = 1 if atom_name in chain_1_atoms else 2 if atom_name in chain_2_atoms else None

        if chain_number:
            output_line = f"{atom_index:5d} {atom_name:5s} {residue_index:5d} {chain_number:5d}\n"
            output_file.write(output_line)

print(f"Data has been written to {output_file_path}")
