import os
import argparse
import sys
import subprocess

def main():
    parser = argparse.ArgumentParser(description="Calculate dihedral angle distributions")
    parser.add_argument("--xtc", help="Path to XTC file", required=True)
    parser.add_argument("--nlip", type=int, default=128, help="Number of lipids (default: 128)")
    parser.add_argument("--offset", type=int, default=134, help="Atom offset per lipid (default: 134 for POPC)")
    parser.add_argument("--out_dir", default="xvg_files", help="Output directory for XVG files", required=True)
    parser.add_argument("--idx_dir", default="index_files", help="Directory for index files", required=True)
    parser.add_argument("--torsion_set", choices=["phosphate", "sn1_isomerization", "double_bond", "all"])
    parser.add_argument("--xvg-suffix", default="_OFF", help="Suffix for XVG output files (default: _OFF)")
    parser.add_argument("--dry-run", action="store_true") #debug
    args = parser.parse_args()

    os.makedirs(args.idx_dir, exist_ok=True)
    os.makedirs(args.out_dir, exist_ok=True)
    
    torsions = get_torsions(args.torsion_set)
    
    process_torsions(torsions, args.nlip, args.offset, args.idx_dir, args.out_dir, args.xtc, args.xvg_suffix, args.dry_run)
    print(f"Done, saved to {args.out_dir}")

def get_torsions(torsion_set):
    
    phosphate_torsions = {
        'P1_O1_C2_H2': [15364, 15363, 15362, 15414], # phosphate torsions
        'P1_O4_C3_H4': [15364, 15367, 15368, 15416],
        'C2_O1_P1_O4': [15362, 15363, 15364, 15367],
        'C2_O1_P1_O3': [15362, 15363, 15364, 15366],
        'C1_C2_O1_P1': [15361, 15362, 15363, 15364],
        'C2_O1_P1_O2': [15362, 15363, 15364, 15365],
        'P1_O1_C2_H3': [15364, 15363, 15362, 15415],
        'P1_O4_C3_C4': [15364, 15367, 15368, 15369],
        'P1_O4_C3_H5': [15364, 15367, 15368, 15417],
    }

    sn1_isomerization = {
        '27_28_29_30': [15397, 15398, 15399, 15400], # sn-1 isomerization torsions
        '28_29_30_31': [15398, 15399, 15400, 15401],
        '29_30_31_32': [15399, 15400, 15401, 15402],
        '30_31_32_33': [15400, 15401, 15402, 15403],
        '31_32_33_34': [15401, 15402, 15403, 15404],
        '32_33_34_35': [15402, 15403, 15404, 15405],
        '33_34_35_36': [15403, 15404, 15405, 15406],
        '34_35_36_37': [15404, 15405, 15406, 15407],
        '35_36_37_38': [15405, 15406, 15407, 15408],
        '36_37_38_39': [15406, 15407, 15408, 15409],
        '37_38_39_40': [15407, 15408, 15409, 15410],
        '38_39_40_41': [15408, 15409, 15410, 15411],
        '39_40_41_42': [15409, 15410, 15411, 15412],
    }
    
    double_bond_torsions = {
        'DB':       [15383, 15384, 15385, 15386], # sn-2 carbon torsions (+adjacent compensatory torsions), nitrogen headgroup, and glycerol backbone
        'DB_HH':    [15443, 15384, 15385, 15444],
        'DB_HC':    [15443, 15384, 15385, 15386],
        '2-5_HH':   [15429, 15377, 15378, 15431],
        '2-5_HC':   [15429, 15377, 15378, 15379],
        '2-5_sn2':  [15377, 15378, 15379, 15380],
        '15-18_sn2':[15390, 15391, 15392, 15393],
        'NHG':      [15370, 15369, 15368, 15367],
        'GB':       [15362, 15361, 15395, 15396],
    }
    
    if torsion_set == "sn1_isomerization":
        return sn1_isomerization
    elif torsion_set == "phosphate":
        return phosphate_torsions
    elif torsion_set == "double_bond":
        return double_bond_torsions
    elif torsion_set == "all":
        all_torsions = {}
        all_torsions.update(sn1_isomerization)
        all_torsions.update(phosphate_torsions)
        all_torsions.update(double_bond_torsions)
        return all_torsions
    else:
        print(f"Unknown torsion set: {torsion_set}")
        sys.exit(1)

def process_torsions(torsions, nlip, offset, idx_dir, out_dir, xtc_file, xvg_suffix="_OFF", dry_run=False):
    successful = 0 # debug, not necessary for read
    failed = 0
    
    for name, atoms in torsions.items():
        idx_file = os.path.join(idx_dir, f"{name}.ndx")
        label = "[ " + name + " ]"

        with open(idx_file, "w") as pfile:
            print(label, file=pfile)
            for n in range(nlip):
                print(f"{atoms[0] + offset * n:4d} {atoms[1] + offset * n:4d} "
                      f"{atoms[2] + offset * n:4d} {atoms[3] + offset * n:4d}", file=pfile)

        xvg_file = os.path.join(out_dir, f"{name}{xvg_suffix}.xvg")
        gmx_command = f"gmx angle -f {xtc_file} -n {idx_file} -od {xvg_file} -type dihedral"
        
        # print(f"Running: {gmx_command}")
        
        if not dry_run:
            try:
                result = subprocess.run(gmx_command, shell=True, check=False)
                
                if result.returncode == 0:
                    print(f"Successfully generated {xvg_file}")
                    successful += 1
                else:
                    print(f"Warning: Command failed with exit code {result.returncode}")
                    failed += 1
            except Exception as e:
                print(f"Error running command: {e}")
                failed += 1
        else:
            print("debug ran")
    
    # print(f"{successful} successful, {failed} failed")

if __name__ == "__main__":
    main()