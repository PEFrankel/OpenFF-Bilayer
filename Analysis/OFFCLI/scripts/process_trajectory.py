import argparse
import os
import sys
import tempfile
import subprocess


def process_trajectory(tpr_file, gro_file=None, xtc_file=None, output_dir="./processed",
                       begin_time=None, skip=None, lipid_selection=None):
#Validation
    if not os.path.exists(tpr_file):
        raise FileNotFoundError(f"TPR file not found: {tpr_file}")
    
    if gro_file is None and xtc_file is None:
        raise ValueError("At least one of gro_file or xtc_file must be provided")
    
    if gro_file and not os.path.exists(gro_file):
        raise FileNotFoundError(f"GRO file not found: {gro_file}")
    
    if xtc_file and not os.path.exists(xtc_file):
        raise FileNotFoundError(f"XTC file not found: {xtc_file}")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Convert paths to absolute paths
    tpr_file = os.path.abspath(tpr_file)
    if gro_file:
        gro_file = os.path.abspath(gro_file)
    if xtc_file:
        xtc_file = os.path.abspath(xtc_file)
    
    # Create temporary selection files that will be used multiple times
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as system_sel:
        system_sel.write("System\n")
        system_sel_path = system_sel.name

    if not lipid_selection:
        lipid_selection = "Lipid"
    
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as lipid_sel:
        lipid_sel.write(f"{lipid_selection}\n")
        lipid_sel.write("System\n")
        lipid_sel_path = lipid_sel.name
    
    # Dictionary to store output file paths
    output_files = {"gro": None, "xtc": None}
    
    # Process GRO file if provided
    if gro_file:
        print(f"\nProcessing GRO file: {gro_file}")
        
        # Step 1: Remove jumps across PBC
        nojump_gro = os.path.join(output_dir, "nojump.gro")
        print(f"Running: gmx trjconv -s {tpr_file} -f {gro_file} -pbc nojump -o {nojump_gro}")
        
        subprocess.run(
            f"gmx trjconv -s {tpr_file} -f {gro_file} -pbc nojump -o {nojump_gro}",
            shell=True, 
            stdin=open(system_sel_path, 'r'),
            check=True
        )
        
        # Step 2: Center on lipids
        centered_gro = os.path.join(output_dir, "centered.gro")
        print(f"Running: gmx trjconv -s {tpr_file} -f {nojump_gro} -pbc mol -center -o {centered_gro}")
        
        subprocess.run(
            f"gmx trjconv -s {tpr_file} -f {nojump_gro} -pbc mol -center -o {centered_gro}",
            shell=True,
            stdin=open(lipid_sel_path, 'r'),
            check=True
        )
        
        # Step 3: Make molecules whole
        analysis_gro = os.path.join(output_dir, "analysis.gro")
        print(f"Running: gmx trjconv -s {tpr_file} -f {centered_gro} -pbc whole -o {analysis_gro}")
        
        subprocess.run(
            f"gmx trjconv -s {tpr_file} -f {centered_gro} -pbc whole -o {analysis_gro}",
            shell=True,
            stdin=open(system_sel_path, 'r'),
            check=True
        )
        
        output_files["gro"] = analysis_gro
        print(f"GRO file processing complete. Final output: {analysis_gro}")
    
    # Process XTC file if provided
    if xtc_file:
        print(f"\nProcessing XTC file: {xtc_file}")
        
        # Build optional parameters
        time_param = f"-b {begin_time}" if begin_time is not None else ""
        skip_param = f"-skip {skip}" if skip is not None else ""
        
        # Step 4: Remove jumps across PBC for trajectory
        nojump_xtc = os.path.join(output_dir, "nojump.xtc")
        cmd = f"gmx trjconv -s {tpr_file} -f {xtc_file} {time_param} {skip_param} -pbc nojump -o {nojump_xtc}"
        print(f"Running: {cmd}")
        
        subprocess.run(
            cmd,
            shell=True,
            stdin=open(system_sel_path, 'r'),
            check=True
        )
        
        # Step 5: Center on lipids
        centered_xtc = os.path.join(output_dir, "centered.xtc")
        print(f"Running: gmx trjconv -s {tpr_file} -f {nojump_xtc} -pbc mol -center -o {centered_xtc}")
        
        subprocess.run(
            f"gmx trjconv -s {tpr_file} -f {nojump_xtc} -pbc mol -center -o {centered_xtc}",
            shell=True,
            stdin=open(lipid_sel_path, 'r'),
            check=True
        )
        
        # Step 6: Make molecules whole
        analysis_xtc = os.path.join(output_dir, "analysis.xtc")
        print(f"Running: gmx trjconv -s {tpr_file} -f {centered_xtc} -pbc whole -o {analysis_xtc}")
        
        subprocess.run(
            f"gmx trjconv -s {tpr_file} -f {centered_xtc} -pbc whole -o {analysis_xtc}",
            shell=True,
            stdin=open(system_sel_path, 'r'),
            check=True
        )
        
        output_files["xtc"] = analysis_xtc
        print(f"XTC file processing complete. Final output: {analysis_xtc}")
    
    os.unlink(system_sel_path)
    os.unlink(lipid_sel_path)
    
    return output_files


def main():
    parser = argparse.ArgumentParser(description="Process GROMACS trajectory files for analysis")
    parser.add_argument("--tpr", required=True, help="Path to TPR topology file")
    parser.add_argument("--gro", help="Path to GRO structure file")
    parser.add_argument("--xtc", help="Path to XTC trajectory file")
    parser.add_argument("--output", default="./data", help="Directory to save processed files")
    parser.add_argument("--begin-time", type=int, help="Beginning time (ps) for trajectory processing")
    parser.add_argument("--skip", type=int, help="Frame skip rate for trajectory processing")
    parser.add_argument("--resname", help="Standard residue name for lipid selection (e.g., POPC)")
    parser.add_argument("--lipid-selection", help="Custom lipid selection string")
    
    args = parser.parse_args()
    
    # Validate inputs
    if not args.gro and not args.xtc:
        parser.error("At least one of --gro or --xtc must be provided")
    
    # Create lipid selection string - prioritize direct selection if provided
    lipid_selection = None
    if args.lipid_selection:
        lipid_selection = args.lipid_selection
    elif args.resname:
        lipid_selection = f"resname {args.resname}"
    
    try:
        output_files = process_trajectory(
            tpr_file=args.tpr,
            gro_file=args.gro,
            xtc_file=args.xtc,
            output_dir=args.output,
            begin_time=args.begin_time,
            skip=args.skip,
            lipid_selection=lipid_selection
        )
        
        print("\nTrajectory processing completed successfully.")
        print(f"Processed GRO: {output_files['gro'] if output_files['gro'] else 'Not processed'}")
        print(f"Processed XTC: {output_files['xtc'] if output_files['xtc'] else 'Not processed'}")
        
    except Exception as e:
        print(f"Error during trajectory processing: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()