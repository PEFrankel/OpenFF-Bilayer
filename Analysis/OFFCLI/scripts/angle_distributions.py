import os
import argparse
import sys
import subprocess

def main():
    parser = argparse.ArgumentParser(description="Calculate angle distributions")
    parser.add_argument("--xtc", help="Path to XTC file", required=True)
    parser.add_argument("--nlip", type=int, default=128, help="Number of lipids (default: 128)")
    parser.add_argument("--offset", type=int, default=134, help="Atom offset per lipid (default: 134 for POPC)")
    parser.add_argument("--out_dir", default="xvg_files", help="Output directory for XVG files")
    parser.add_argument("--idx_dir", default="index_files", help="Directory for index files")
    parser.add_argument("--angle_set", choices=["phosphate", "standard", "all"])
    parser.add_argument("--xvg-suffix", default="_OFF", help="Suffix for XVG output files (default: _OFF)")
    parser.add_argument("--dry-run", action="store_true") #debug
    args = parser.parse_args()

    os.makedirs(args.idx_dir, exist_ok=True)
    os.makedirs(args.out_dir, exist_ok=True)
    
    angles = get_angles(args.angle_set)

    process_angles(angles, args.nlip, args.offset, args.idx_dir, args.out_dir, args.xtc, args.xvg_suffix, args.dry_run)
    print(f"Done, saved to {args.out_dir}")

def get_angles(angle_set="phosphate"):
    phosphate_angles = {
        'O1_P1_O4': [15363, 15364, 15367], # phosphate angles (add SMIRKS IDs, only t159 and t160 here)
        'O2_P1_O3': [15365, 15364, 15366],
        'O3_P1_O4': [15366, 15364, 15367],
        'O1_P1_O3': [15363, 15364, 15366],
        'O1_P1_O2': [15363, 15364, 15365],
        'C2_O1_P1': [15362, 15363, 15364],
        'O2_P1_O4': [15365, 15364, 15367],
    }
    
    standard_angles = {
        'C1_C2_C3': [15361, 15362, 15368],  # glycerol backbone angles
        'C2_C3_O4': [15362, 15368, 15367],
        'C3_C2_O1': [15368, 15362, 15363],
        'N1_C4_C5': [15370, 15371, 15372],  # Amide headgroup angles
    }
    
    if angle_set == "phosphate":
        return phosphate_angles
    elif angle_set == "standard":
        return standard_angles
    elif angle_set == "all":
        all_angles = {}
        all_angles.update(phosphate_angles)
        all_angles.update(standard_angles)
        return all_angles
    else:
        print(f"Unknown angle set: {angle_set}")
        sys.exit(1)

def process_angles(angles, nlip, offset, idx_dir, out_dir, xtc_file, xvg_suffix="_OFF", dry_run=False):
    successful = 0 # debug, not necessary for read
    failed = 0
    
    for name, atoms in angles.items():
        idx_file = os.path.join(idx_dir, f"{name}.ndx")
        label = "[ " + name + " ]"

        with open(idx_file, "w") as pfile:
            print(label, file=pfile)
            for n in range(nlip):
                print(f"{atoms[0] + offset * n:4d} {atoms[1] + offset * n:4d} {atoms[2] + offset * n:4d}", file=pfile)

        xvg_file = os.path.join(out_dir, f"{name}{xvg_suffix}.xvg")
        gmx_command = f"gmx angle -f {xtc_file} -n {idx_file} -od {xvg_file} -type angle"
        
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
    
    print(f"{successful} successful, {failed} failed")

if __name__ == "__main__":
    main()