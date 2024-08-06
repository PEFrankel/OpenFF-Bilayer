from rdkit import Chem
from rdkit.Chem import AllChem

# cholesterol: C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C
# ergosterol: 'C[C@H](/C=C/[C@H](C)C(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3C2=CC=C4[C@@]3(CC[C@@H](C4)O)C)C
# beta-Sitosterol: CC[C@H](CC[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C)C(C)C
# stigmasterol : CC[C@H](/C=C/[C@@H](C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C)C(C)C
# cortisol: [3H]C1[C@H]2[C@@H]3CC[C@@]([C@]3(C[C@@H]([C@@H]2[C@]4(C(C(C(=O)C=C4C1[3H])[3H])[3H])C)O)C)(C(=O)CO)O
# corticosterone: C[C@]12CCC(=O)C=C1CC[C@@H]3[C@@H]2[C@H](C[C@]4([C@H]3CC[C@@H]4C(=O)CO)C)O
# aldosterone: C[C@]12CCC(=O)C=C1CC[C@@H]3[C@@H]2[C@H](C[C@]4([C@H]3CC[C@@H]4C(=O)CO)C=O)O
# progesterone: CC(=O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4=CC(=O)CC[C@]34C)C
# beta-estradiol: C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O
# testosterone: C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@]34C
# estradiol: C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=C3C=CC(=C4)O

def parse_smiles(smiles):

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    return mol

def remove_fragment(mol, fragment_smiles):

    print(f"Removing fragment: {fragment_smiles}")
    fragment = Chem.MolFromSmiles(fragment_smiles)
    if fragment is None:
        raise ValueError(f"Invalid SMILES string: {fragment_smiles}")
    
    matches = mol.GetSubstructMatches(fragment)
    if not matches:
        raise ValueError(f"Fragment not found: {fragment_smiles}")
    
    print(f"Found fragment matches: {matches}")

    # new molecule without the matched fragment
    editable_mol = Chem.EditableMol(mol)
    atoms_to_remove = set([atom_idx for match in matches for atom_idx in match])
    for atom_idx in sorted(atoms_to_remove, reverse=True):
        editable_mol.RemoveAtom(atom_idx)
    
    new_mol = editable_mol.GetMol()
    Chem.SanitizeMol(new_mol)
    print(f"molecule after fragment removal: {Chem.MolToSmiles(new_mol)}")
    return new_mol

def get_attachment_point(mol, substructure_smiles):

    print(f"Attachment point for: {substructure_smiles}")
    pattern = Chem.MolFromSmiles(substructure_smiles)
    if pattern is None:
        raise ValueError(f"Invalid SMILES string: {substructure_smiles}")
    
    matches = mol.GetSubstructMatches(pattern)
    if matches:
        print(f"Found matches for {substructure_smiles}: {matches}")
        return matches[0][0]
    raise ValueError(f"No attachment found for {substructure_smiles}")

def combine_molecules(tail_mol, steroid_mol, tail_attach_idx, steroid_attach_idx):

    combined_mol = Chem.RWMol(Chem.CombineMols(tail_mol, steroid_mol))
    steroid_attach_idx += tail_mol.GetNumAtoms()
    
    # need to add a carbon/bond at attachment for just this cort/chol connection
    combined_mol.AddBond(tail_attach_idx, steroid_attach_idx, Chem.BondType.SINGLE)
    
    # need to remove hydrogen for just this case again
    atom = combined_mol.GetAtomWithIdx(steroid_attach_idx)
    for neighbor in atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 1:  # remove hydrogen atoms only
            combined_mol.RemoveAtom(neighbor.GetIdx())
            break
    
    new_mol = combined_mol.GetMol()
    Chem.SanitizeMol(new_mol)
    return new_mol

def validate_smiles(smiles):

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Invalid SMILES string: {smiles}")
    return mol is not None

def main():
    try:
        # test with chol and cort
        cholesterol = "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C"
        corticosterone = "C[C@]12CCC(=O)C=C1CC[C@@H]3[C@@H]2[C@H](C[C@]4([C@H]3CC[C@@H]4)C)O"
        
        # tranlsate SMILES strings into RDKit objects
        cholesterol_mol = parse_smiles(cholesterol)
        corticosterone_mol = parse_smiles(corticosterone)
        
        # separate tail from steroid
        cholesterol_tail_smiles = "C[C@H](CCCC(C)C)"
        corticosterone_steroid_smiles = "C[C@]12CCC(=O)C=C1CC[C@@H]3[C@@H]2[C@H](C[C@]4([C@H]3CC[C@@H]4)C)O"
        
        cholesterol_tail = Chem.MolFromSmiles(cholesterol_tail_smiles)
        print(f"Cholesterol tail: {Chem.MolToSmiles(cholesterol_tail)}")
        corticosterone_steroid = remove_fragment(corticosterone_mol, cholesterol_tail_smiles)
        
        print(f"Corticosterone steroid after removal: {Chem.MolToSmiles(corticosterone_steroid)}")
        
        tail_attach_idx = get_attachment_point(cholesterol_mol, cholesterol_tail_smiles)
        steroid_attach_idx = get_attachment_point(corticosterone_steroid, corticosterone_steroid_smiles)
        
        # combine
        new_mol = combine_molecules(cholesterol_tail, corticosterone_steroid, tail_attach_idx, steroid_attach_idx)
        
        # read new smiles string
        new_smiles = Chem.MolToSmiles(new_mol, isomericSmiles=True)
        
        # validate the smiles string
        is_valid = validate_smiles(new_smiles)
        
        if is_valid:
            print("Valid SMILES string:", new_smiles)
        else:
            print("No valid SMILES string.")
    
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
