import os

nlip = 128
offset = 134

#  change indices/torsion names here if needed
torsions = {
    ## OpenFF (see original 2.2.0 SMILES string)
    'O1_P1_O4': [15363, 15364, 15367], # phosphate angles (SMIRKS IDs: t159 and t160)
    'O2_P1_O3': [15365, 15364, 15366],
    'O3_P1_O4': [15366, 15364, 15367],
    'O1_P1_O3': [15363, 15364, 15366],
    'O1_P1_O2': [15363, 15364, 15365],
    'C2_O1_P1': [15362, 15363, 15364],
    'O2_P1_O4': [15365, 15364, 15367],
    # 'C1_C2_C3': [15361, 15362, 15368],  # glycerol backbone angles
    # 'C2_C3_O4': [15362, 15368, 15367],
    # 'C3_C2_O1': [15368, 15362, 15363],
    # 'N1_C4_C5': [15370, 15371, 15372],  # Amide headgroup angles
    ## Lipid21
    # 'O1_P1_O4': [58, 59, 60],
    # 'O2_P1_O3': [81, 59, 80],
    # 'O3_P1_O4': [80, 59, 60],
    # 'O1_P1_O3': [58, 59, 80],
    # 'O1_P1_O2': [58, 59, 81],
    # 'C2_O1_P1': [55, 58, 59],
    # 'O2_P1_O4': [81, 59, 60],
    ## CHARMM36
    # 'O1_P1_O4': [24, 20, 23],
    # 'O2_P1_O3': [22, 20, 21],
    # 'O3_P1_O4': [21, 20, 23],
    # 'O1_P1_O3': [24, 20, 21],
    # 'O1_P1_O2': [24, 20, 22],
    # 'C2_O1_P1': [25, 24, 20],
    # 'O2_P1_O4': [22, 20, 23],
}

os.makedirs('index_files', exist_ok=True)
os.makedirs('xvg_files', exist_ok=True)

for name, atoms in torsions.items():
    fname = os.path.join('index_files', name + ".ndx")
    label = "[ " + name + " ]"

    with open(fname, "w") as pfile:
        print(label, file=pfile)
        for n in range(nlip):
            print(f"{atoms[0] + offset * n:4d} {atoms[1] + offset * n:4d} {atoms[2] + offset * n:4d}", file=pfile)

# "_OFF" is the default xvg suffix I have here. Change if needed
    xvg_name = f"xvg_files/{name}_OFF.xvg"
    gmx_command = f"gmx angle -f analysis.xtc -n {fname} -od {xvg_name} -type angle"
    os.system(gmx_command)
