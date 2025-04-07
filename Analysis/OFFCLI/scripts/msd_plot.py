import argparse
import os
import subprocess
import sys
import time

def create_phosphate_index(tpr_path, index_dir, phosphate_atom='P1x'):
    os.makedirs(index_dir, exist_ok=True)
    index_path = os.path.join(index_dir, "phosphates.ndx")

    try:
        subprocess.run(
            f"gmx make_ndx -f {tpr_path} -o {index_path}",
            shell=True,
            input=f'a {phosphate_atom}\nq\n', 
            text=True, 
            check=True
        )
    except subprocess.CalledProcessError as e:
        # print(f"debug, creating phosphate index")
        sys.exit(1)
    
    return index_path

def run_msd_analysis(tpr_path, xtc_path, index_path, xvg_dir, phosphate_atom='P1x', prefix=''):
    os.makedirs(xvg_dir, exist_ok=True)
    output_filename = f"{prefix}msd.xvg" if prefix else "msd.xvg"
    output_path = os.path.join(xvg_dir, output_filename)
    
    print("="*102)
    print("WARNING: Ensure trajectory was processed with -pbc nojump to prevent overestimation of displacement")
    print("="*102 + "\n")
    
    print("-"*102)
    print("5 second limbo for you to recall if PBC was(were) accounted for in case `process_trajectory.py` wasn't used")
    print("-"*102 + "\n"*2)
    time.sleep(8)
    
    try:
        process = subprocess.Popen( # need to be very specific with gmx inputs after intial call
            f"gmx msd -s {tpr_path} -f {xtc_path} -n {index_path} -lateral z -o {output_path}",
            shell=True,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )

        stdout, stderr = process.communicate(input=f'{phosphate_atom}\n')
        
        if process.returncode != 0:
            print(f"error: {stderr}")
            sys.exit(1)
        print(f"Output saved to {output_path}")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Mean Square Displacement of Phosphate Atoms")
    parser.add_argument("--tpr", required=True, help="Path to TPR file")
    parser.add_argument("--xtc", help="Defaults to `analysis.xtc` from process_trajectory.py")
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--index-dir")
    parser.add_argument("--xvg-dir")
    parser.add_argument("--prefix", default="", help="Prefix for output file name")
    parser.add_argument("--phosphate-atom", default='P1x', help="Atom name to use for phosphate selection in index creation (default: 'P1x')")
    
    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    index_dir = args.index_dir or os.path.join(args.out_dir, "index_files")
    xvg_dir = args.xvg_dir or os.path.join(args.out_dir, "xvg_files")
    
    if not args.xtc:
        default_processed_xtc = os.path.join(os.path.dirname(args.out_dir), "process_output", "data", "analysis.xtc")
        if os.path.exists(default_processed_xtc):
            args.xtc = default_processed_xtc
        else:
            print("No trj")
            sys.exit(1)
    
    phosphate_index = create_phosphate_index(args.tpr, index_dir, args.phosphate_atom)
    run_msd_analysis(args.tpr, args.xtc, phosphate_index, xvg_dir, args.phosphate_atom, args.prefix)

if __name__ == "__main__":
    main()