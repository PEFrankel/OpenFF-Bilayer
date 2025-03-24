import argparse
import os
import subprocess
import sys

def main():
    parser = argparse.ArgumentParser(description="Lipid bilayer analysis tool")
    parser.add_argument("--analysis", choices=["process", "dihedral", "angle", "apl", "rdf", "isomerization"], 
                        required=True, help="Type of analysis to perform")
    parser.add_argument("--xtc", help="Path to XTC trajectory file")
    parser.add_argument("--tpr", help="Path to TPR topology file")
    parser.add_argument("--gro", help="Path to GRO structure file")
    parser.add_argument("--nlip", type=int, help="Number of lipids")
    parser.add_argument("--offset", type=int, help="Atom offset per lipid")
    parser.add_argument("--resname", help="Residue name for standard residue analysis")
    parser.add_argument("--lipid21-resnames", nargs='+', help="Residue names for split lipids (e.g., PA PC OL for Lipid17/21)")
    parser.add_argument("--torsion_set", choices=["sn1_isomerization", "phosphate", "double_bond", "all"], 
                        help="Torsion set to analyze for dihedral analysis")
    parser.add_argument("--title", help="Title for the plot and output files")
    parser.add_argument("--max-time", type=int, help="Maximum time (ps) to include in analysis")
    parser.add_argument("--color", help="Color to use for plotting (force field or basic color name)")
    parser.add_argument("--output-dir", default="output", help="Path to output directory")
    parser.add_argument("--scripts-dir", help="Path to directory containing analysis scripts")
    parser.add_argument("--xvg-suffix", default="_aux", help="Suffix to use for XVG output files (default: _aux)")
    
    # Arguments for trajectory processing (from process_trajectory.py)
    parser.add_argument("--begin-time", type=int, help="Beginning time (ps) for trajectory processing")
    parser.add_argument("--skip", type=int, help="Frame skip rate for trajectory processing")
    parser.add_argument("--lipid-selection", help="Selection string for lipids (default: 'Lipid' or based on resname)")

    args = parser.parse_args()

    # Determine scripts directory - either from argument or relative to this script
    if args.scripts_dir:
        scripts_dir = os.path.abspath(args.scripts_dir)
    else:
        # Default to the directory where this script is located
        scripts_dir = os.path.dirname(os.path.abspath(__file__))

    # Create output directory if it doesn't exist
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)
    
    # Create subdirectories for organized output
    xvg_dir = os.path.join(output_dir, "xvg_files")
    idx_dir = os.path.join(output_dir, "index_files")
    os.makedirs(xvg_dir, exist_ok=True)
    os.makedirs(idx_dir, exist_ok=True)

    # Run the appropriate script based on the analysis type
    if args.analysis == "process":
        run_trajectory_processing(args, scripts_dir, output_dir)
    elif args.analysis == "dihedral":
        run_dihedral_analysis(args, scripts_dir, output_dir, xvg_dir, idx_dir)
    elif args.analysis == "angle":
        run_angle_analysis(args, scripts_dir, output_dir, xvg_dir, idx_dir)
    elif args.analysis == "apl":
        run_apl_analysis(args, scripts_dir, output_dir)
    elif args.analysis == "rdf":
        run_rdf_analysis(args, scripts_dir, output_dir)
    elif args.analysis == "isomerization":
        run_isomerization_analysis(args, scripts_dir, output_dir)

def run_trajectory_processing(args, scripts_dir, output_dir):
    script_path = os.path.join(scripts_dir, "process_trajectory.py")
    cmd = [sys.executable, script_path]
    
    # Check for required arguments
    if not args.tpr:
        print("Error: --tpr is required for trajectory processing")
        sys.exit(1)
    
    cmd.extend(["--tpr", os.path.abspath(args.tpr)])
    
    # At least one of --gro or --xtc must be provided
    if not (args.gro or args.xtc):
        print("Error: At least one of --gro or --xtc must be provided for trajectory processing")
        sys.exit(1)
    
    if args.gro:
        cmd.extend(["--gro", os.path.abspath(args.gro)])
    if args.xtc:
        cmd.extend(["--xtc", os.path.abspath(args.xtc)])
    
    # Process trajectory-specific arguments
    if args.begin_time is not None:
        cmd.extend(["--begin-time", str(args.begin_time)])
    if args.skip is not None:
        cmd.extend(["--skip", str(args.skip)])
    
    # Handle lipid selection - either from explicit lipid-selection or from resname
    if args.lipid_selection:
        cmd.extend(["--lipid-selection", args.lipid_selection])
    elif args.resname:
        cmd.extend(["--resname", args.resname])
    
    # Change this line to direct output to the data directory
    data_dir = os.path.join(output_dir, "data")
    os.makedirs(data_dir, exist_ok=True)
    cmd.extend(["--output", data_dir])
    
    print(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd)

def run_dihedral_analysis(args, scripts_dir, output_dir, xvg_dir, idx_dir):
    script_path = os.path.join(scripts_dir, "dihedral_distributions.py")
    cmd = [sys.executable, script_path]
    
    if not args.xtc:
        print("Error: --xtc is required for dihedral analysis")
        sys.exit(1)
    
    cmd.extend(["--xtc", os.path.abspath(args.xtc)])
    cmd.extend(["--out_dir", xvg_dir])
    cmd.extend(["--idx_dir", idx_dir])
    
    if args.nlip:
        cmd.extend(["--nlip", str(args.nlip)])
    if args.offset:
        cmd.extend(["--offset", str(args.offset)])
    if args.torsion_set:
        cmd.extend(["--torsion_set", args.torsion_set])
    if args.xvg_suffix:
        cmd.extend(["--xvg-suffix", args.xvg_suffix])
    
    print(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd)

def run_angle_analysis(args, scripts_dir, output_dir, xvg_dir, idx_dir):
    script_path = os.path.join(scripts_dir, "angle_distributions.py")
    cmd = [sys.executable, script_path]
    
    if not args.xtc:
        print("Error: --xtc is required for angle analysis")
        sys.exit(1)
    
    cmd.extend(["--xtc", os.path.abspath(args.xtc)])
    cmd.extend(["--out_dir", xvg_dir])
    cmd.extend(["--idx_dir", idx_dir])
    
    if args.nlip:
        cmd.extend(["--nlip", str(args.nlip)])
    if args.offset:
        cmd.extend(["--offset", str(args.offset)])
    if args.xvg_suffix:
        cmd.extend(["--xvg-suffix", args.xvg_suffix])
    
    print(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd)

def run_apl_analysis(args, scripts_dir, output_dir):
    script_path = os.path.join(scripts_dir, "APLfromTRJ.py")
    cmd = [sys.executable, script_path]
    
    if not args.gro:
        print("Error: --gro is required for APL analysis")
        sys.exit(1)
    cmd.extend(["--gro", os.path.abspath(args.gro)])
    if not args.xtc:
        print("Error: --xtc is required for APL analysis")
        sys.exit(1)
    cmd.extend(["--xtc", os.path.abspath(args.xtc)])
    
    # Either standard resname or lipid21 resnames are necessary
    if not (args.resname or args.lipid21_resnames):
        print("Error: Either --resname or --lipid21-resnames must be provided for APL analysis")
        sys.exit(1)
    if args.resname:
        cmd.extend(["--standard-resname", args.resname])
    if args.lipid21_resnames:
        cmd.extend(["--lipid21-resnames"] + args.lipid21_resnames)
    if args.title:
        cmd.extend(["--title", args.title])
    if args.max_time:
        cmd.extend(["--max-time", str(args.max_time)])
    if args.color:
        cmd.extend(["--color", args.color])
        
    # Add output directory
    output_prefix = os.path.join(output_dir, "apl_results")
    cmd.extend(["--output", output_prefix])
    
    print(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd)

def run_rdf_analysis(args, scripts_dir, output_dir):
    script_path = os.path.join(scripts_dir, "RDF_multiprocess_PO-NC_IDexclusion.py")
    cmd = [sys.executable, script_path]
    
    if not args.xtc:
        print("Error: --xtc is required for RDF analysis")
        sys.exit(1)
    cmd.extend(["--xtc", os.path.abspath(args.xtc)])
    
    if not args.tpr:
        print("Error: --tpr is required for RDF analysis")
        sys.exit(1)
    cmd.extend(["--tpr", os.path.abspath(args.tpr)])
    
    # Add output directory and filename
    output_path = os.path.join(output_dir, "rdf_analysis")
    cmd.extend(["--output", output_path])
    
    # Add output directory path for index and xvg files
    index_dir = os.path.join(output_dir, "index_files")
    xvg_dir = os.path.join(output_dir, "xvg_files")
    cmd.extend(["--index-dir", index_dir])
    cmd.extend(["--xvg-dir", xvg_dir])
    
    print(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd)

def run_isomerization_analysis(args, scripts_dir, output_dir):
    script_path = os.path.join(scripts_dir, "isomerization_code_PointGraph.py")
    cmd = [sys.executable, script_path]
    
    # Create data directory for organizing input XVG folders
    data_dir = os.path.join(output_dir, "data")
    os.makedirs(data_dir, exist_ok=True)
    
    # Define the three XVG folders within the data directory
    xvg_folders = [
        os.path.join(data_dir, "xvg_files_OFF_HMR"),
        os.path.join(data_dir, "xvg_files_macrog"), 
        os.path.join(data_dir, "xvg_files_slipids")
    ]
    
    # Create the folders
    for folder in xvg_folders:
        os.makedirs(folder, exist_ok=True)
        
    # Set the output file path for the plot
    output_file = os.path.join(output_dir, "sn1_isomerization_comparison.png")
    
    # Use the correct parameter names as expected by isomerization_code_PointGraph.py
    cmd.extend(["--xvg_folders"] + xvg_folders)
    cmd.extend(["--output", output_file])
    
    if args.title:
        cmd.extend(["--title", args.title])
    
    print(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd)

if __name__ == "__main__":
    main()