import os
import argparse
import sys
import subprocess

def main():
    parser = argparse.ArgumentParser(description="Calculate angle distributions from MD trajectories")
    parser.add_argument("--xtc", help="Trajectory file (.xtc format)", required=True)
    parser.add_argument("--nlip", type=int, default=128, help="Number of lipids (default: 128)")
    parser.add_argument("--offset", type=int, default=134, help="Atom offset per lipid (default: 134 for POPC)")
    parser.add_argument("--out_dir", default="xvg_files", help="Output directory for XVG files")
    parser.add_argument("--idx_dir", default="index_files", help="Directory for index files")
    parser.add_argument("--angle_set", choices=["phosphate", "standard", "all"], 
                        default="phosphate", help="Angle set to analyze (default: phosphate)")
    parser.add_argument("--xvg-suffix", default="_aux", help="Suffix to use for XVG output files (default: _aux)")
    parser.add_argument("--dry-run", action="store_true", help="Print commands without executing them")
    
    args = parser.parse_args()
    
    # Check if trajectory file exists
    if not os.path.exists(args.xtc):
        print(f"Error: Trajectory file '{args.xtc}' does not exist")
        sys.exit(1)
    
    # Create output directories
    os.makedirs(args.idx_dir, exist_ok=True)
    os.makedirs(args.out_dir, exist_ok=True)
    
    # Get angle definitions based on the requested set
    angles = get_angles(args.angle_set)
    
    # Process the angles
    process_angles(angles, args.nlip, args.offset, args.idx_dir, args.out_dir, args.xtc, args.xvg_suffix, args.dry_run)
    
    print(f"Angle distribution analysis complete. Results saved in {args.out_dir}")

def get_angles(angle_set="phosphate"):
    """Get the angle definitions based on the requested set"""
    phosphate_angles = {
        'O1_P1_O4': [15363, 15364, 15367], # Phosphate Angles (add SMIRKS IDs, only t159 and t160 here)
        'O2_P1_O3': [15365, 15364, 15366],
        'O3_P1_O4': [15366, 15364, 15367],
        'O1_P1_O3': [15363, 15364, 15366],
        'O1_P1_O2': [15363, 15364, 15365],
        'C2_O1_P1': [15362, 15363, 15364],
        'O2_P1_O4': [15365, 15364, 15367],
    }
    
    standard_angles = {
        'C1_C2_C3': [15361, 15362, 15368],  # Glycerol backbone angles
        'C2_C3_O4': [15362, 15368, 15367],
        'C3_C2_O1': [15368, 15362, 15363],
        'N1_C4_C5': [15370, 15371, 15372],  # Headgroup angles
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

def process_angles(angles, nlip, offset, idx_dir, out_dir, xtc_file, xvg_suffix="_aux", dry_run=False):
    """Process each angle by creating index files and running gmx angle"""
    print(f"Processing {len(angles)} angles with {nlip} lipids (offset: {offset})")
    
    successful = 0
    failed = 0
    
    for name, atoms in angles.items():
        # Create index files
        idx_file = os.path.join(idx_dir, f"{name}.ndx")
        label = "[ " + name + " ]"

        with open(idx_file, "w") as pfile:
            print(label, file=pfile)
            for n in range(nlip):
                print(f"{atoms[0] + offset * n:4d} {atoms[1] + offset * n:4d} {atoms[2] + offset * n:4d}", file=pfile)

        # Run gmx angle
        xvg_file = os.path.join(out_dir, f"{name}{xvg_suffix}.xvg")
        gmx_command = f"gmx angle -f {xtc_file} -n {idx_file} -od {xvg_file} -type angle"
        
        print(f"Running: {gmx_command}")
        
        if not dry_run:
            try:
                # Using subprocess for better error handling
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
            print("Dry run mode: command not executed")
    
    print(f"Angle analysis completed: {successful} successful, {failed} failed")

if __name__ == "__main__":
    main()