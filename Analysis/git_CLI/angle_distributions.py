import os
import argparse
import sys

def main():
    parser = argparse.ArgumentParser(description="Calculate angle distributions from MD trajectories")
    parser.add_argument("--xtc", help="Trajectory file (.xtc format)", required=True)
    parser.add_argument("--nlip", type=int, default=128, help="Number of lipids (default: 128)")
    parser.add_argument("--offset", type=int, default=134, help="Atom offset per lipid (default: 134 for POPC)")
    parser.add_argument("--out_dir", default="xvg_files", help="Output directory for XVG files (default: xvg_files)")
    parser.add_argument("--idx_dir", default="index_files", help="Directory for index files (default: index_files)")
    
    args = parser.parse_args()
    
    if "TRAJECTORY_FILE" in os.environ and not args.xtc:
        args.xtc = os.environ["TRAJECTORY_FILE"]
    if "NLIP" in os.environ and not args.nlip:
        args.nlip = int(os.environ["NLIP"])
    if "OFFSET" in os.environ and not args.offset:
        args.offset = int(os.environ["OFFSET"])
    
    os.makedirs(args.idx_dir, exist_ok=True)
    os.makedirs(args.out_dir, exist_ok=True)
    
    # default OpenFF SMILES
    angles = get_angles()
    
    process_angles(angles, args.nlip, args.offset, args.idx_dir, args.out_dir, args.xtc)

def get_angles():
    """Get the phosphate angle definitions using the default OpenFF SMILES"""
    return {
        'O1_P1_O4': [15363, 15364, 15367], # Phosphate Angles (add SMIRKS IDs, only t159 and t160 here)
        'O2_P1_O3': [15365, 15364, 15366],
        'O3_P1_O4': [15366, 15364, 15367],
        'O1_P1_O3': [15363, 15364, 15366],
        'O1_P1_O2': [15363, 15364, 15365],
        'C2_O1_P1': [15362, 15363, 15364],
        'O2_P1_O4': [15365, 15364, 15367],
    }

def process_angles(angles, nlip, offset, idx_dir, out_dir, xtc_file):
    """Process each angle by creating index files and running gmx angle"""
    print(f"Processing {len(angles)} angles with {nlip} lipids (offset: {offset})")
    
    for name, atoms in angles.items():
        # Create index files
        fname = os.path.join(idx_dir, name + ".ndx")
        label = "[ " + name + " ]"

        with open(fname, "w") as pfile:
            print(label, file=pfile)
            for n in range(nlip):
                print(f"{atoms[0] + offset * n:4d} {atoms[1] + offset * n:4d} {atoms[2] + offset * n:4d}", file=pfile)

        # Run gmx angle
        xvg_name = os.path.join(out_dir, f"{name}_aux.xvg")
        gmx_command = f"gmx angle -f {xtc_file} -n {fname} -od {xvg_name} -type angle"
        
        print(f"Running: {gmx_command}")
        result = os.system(gmx_command)
        
        if result != 0:
            print(f"Warning: Command failed with exit code {result}")
        else:
            print(f"Successfully generated {xvg_name}")

if __name__ == "__main__":
    main()