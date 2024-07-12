import os

nlip = 128
offset = 134

#  change indices/torsion names here if needed
torsions = {
    'P1_O1_C2_H2': [15364, 15363, 15362, 15414],
    'P1_O4_C3_H4': [15364, 15367, 15368, 15416],
    'C2_O1_P1_O4': [15362, 15363, 15364, 15367],
    'C2_O1_P1_O3': [15362, 15363, 15364, 15366],
    'C1_C2_O1_P1': [15361, 15362, 15363, 15364],
    'C2_O1_P1_O2': [15362, 15363, 15364, 15365],
    'P1_O1_C2_H3': [15364, 15363, 15362, 15415],
    'P1_O4_C3_C4': [15364, 15367, 15368, 15369],
    'P1_O4_C3_H5': [15364, 15367, 15368, 15417],
}

os.makedirs('index_files', exist_ok=True)
os.makedirs('xvg_files', exist_ok=True)

for name, atoms in torsions.items():
    fname = os.path.join('index_files', name + ".ndx")
    label = "[ " + name + " ]"

    with open(fname, "w") as pfile:
        print(label, file=pfile)
        for n in range(nlip):
            print(f"{atoms[0] + offset * n:4d} {atoms[1] + offset * n:4d} {atoms[2] + offset * n:4d} {atoms[3] + offset * n:4d}", file=pfile)

    xvg_name = f"xvg_files/{name}_OFF_auto.xvg"
    gmx_command = f"gmx angle -f md.xtc -n {fname} -od {xvg_name} -type dihedral"
    os.system(gmx_command)
