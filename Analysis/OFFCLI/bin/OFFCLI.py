import argparse
import os
import subprocess
import sys

def main():
    parser = argparse.ArgumentParser(description="Lipid bilayer analysis tool")
    parser.add_argument("--analysis", choices=["process", "dihedral", "angle", "apl", "rdf", "isomerization", "msd", "diffusion"], 
                        required=True, help="Type of analysis to perform")
    
    parser.add_argument("--xtc", help="Path to XTC file")
    parser.add_argument("--tpr", help="Path to TPR file")
    parser.add_argument("--gro", help="Path to GRO file")
    parser.add_argument("--xvg", help="Path to XVG file")
    #apl
    parser.add_argument("--nlip", type=int, help="Number of lipids")
    parser.add_argument("--offset", type=int, help="Atom offset per lipid in the gro file")
    parser.add_argument("--resname")
    parser.add_argument("--lipid21-resnames", nargs='+', help="Residue names for split lipids (e.g. PA PC OL for Lipid17/21)")
    parser.add_argument("--max-time", type=int, help="Maximum time (ps) to include in analysis. You probably don't need this")
    # auto torsions/angles
    parser.add_argument("--torsion_set", choices=["sn1_isomerization", "phosphate", "double_bond", "all"])
    parser.add_argument("--angle_set", choices=["phosphate", "standard", "all"])
    parser.add_argument("--title", help="Title for the plot and output files")
    parser.add_argument("--color", help="Color to use for plotting")
    parser.add_argument("--base-dir", default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))), help="Base directory for the project (default: parent of script location)")
    parser.add_argument("--xvg-suffix", default="_OFF", help="Suffix for XVG output files (default: _OFF)")
    # xtc/gro processing
    parser.add_argument("--begin-time", type=int, help="Start time (ps) for trj processing")
    parser.add_argument("--skip", type=int, help="Frame skip for trajectory processing")
    parser.add_argument("--lipid-selection", help="Alt election for lipids (default: 'Lipid' or resname)")
    # msd
    parser.add_argument("--prefix", help="Prefix for MSD output file")
    parser.add_argument("--phosphate-atom", default='P1x', help="Atom name to use for phosphate selection (default: 'P1x')")     
    # diffusion
    parser.add_argument("--window-size", type=int, default=10000, help="Default: 10000 = 100 ns")
    parser.add_argument("--window-step", type=int, default=200, help="Default: 200 = 2 ns")
    parser.add_argument("--min-r2", type=float, default=0.99, help="Minimum linear fit value (default: 0.99)")
    args = parser.parse_args()

    base_dir = os.path.abspath(args.base_dir)
    inputs_dir = os.path.join(base_dir, "inputs")
    outputs_dir = os.path.join(base_dir, "outputs")
    scripts_dir = os.path.join(base_dir, "scripts")
    out_dir = os.path.join(outputs_dir, f"{args.analysis}_output")
    os.makedirs(out_dir, exist_ok=True)

    def resolve_input_path(input_type, filename):
        return os.path.join(inputs_dir, input_type, filename) if filename else None
    
     # checks processed ver first
    xtc_path = None
    if args.xtc:
        processed_xtc = os.path.join(outputs_dir, "process_output", "data", "analysis.xtc")# (need to rmv redundant data dir - artifact)
        if os.path.exists(processed_xtc):
            xtc_path = processed_xtc
        else:
            xtc_path = resolve_input_path("trajectories", os.path.basename(args.xtc))

    tpr_path = resolve_input_path("topologies", os.path.basename(args.tpr)) if args.tpr else None
    gro_path = resolve_input_path("structures", os.path.basename(args.gro)) if args.gro else None # the processed ver of this doesn't change data or save time
    
    xvg_path = None
    if args.xvg:
        if os.path.exists(args.xvg):
            xvg_path = args.xvg
        else:
            xvg_path = resolve_input_path("xvg_files", os.path.basename(args.xvg))

    if args.analysis == "dihedral":
        xvg_dir = os.path.join(out_dir, "xvg_files")
        idx_dir = os.path.join(out_dir, "index_files")
        os.makedirs(xvg_dir, exist_ok=True)
        os.makedirs(idx_dir, exist_ok=True)
        run_dihedral_analysis(args, scripts_dir, out_dir, xvg_dir, idx_dir, xtc_path)
    elif args.analysis == "angle":
        xvg_dir = os.path.join(out_dir, "xvg_files")
        idx_dir = os.path.join(out_dir, "index_files")
        os.makedirs(xvg_dir, exist_ok=True)
        os.makedirs(idx_dir, exist_ok=True)
        run_angle_analysis(args, scripts_dir, out_dir, xvg_dir, idx_dir, xtc_path)
    elif args.analysis == "apl":
        plot_dir = os.path.join(out_dir, "plot_outputs")
        os.makedirs(plot_dir, exist_ok=True)
        run_apl_analysis(args, scripts_dir, out_dir, gro_path, xtc_path)
    elif args.analysis == "rdf":
        index_dir = os.path.join(out_dir, "index_files")
        xvg_dir = os.path.join(out_dir, "xvg_files")
        os.makedirs(index_dir, exist_ok=True)
        os.makedirs(xvg_dir, exist_ok=True)
        run_rdf_analysis(args, scripts_dir, out_dir, xtc_path, tpr_path)
    elif args.analysis == "isomerization":
        data_dir = os.path.join(out_dir, "data")
        plot_dir = os.path.join(out_dir, "plot_outputs")
        os.makedirs(data_dir, exist_ok=True)
        os.makedirs(plot_dir, exist_ok=True)
        xvg_dirs = [                                     # comparison ffs to be updated
            os.path.join(data_dir, "xvg_files_OFF_HMR"),
            os.path.join(data_dir, "xvg_files_macrog"), 
            os.path.join(data_dir, "xvg_files_slipids")
        ]
        
        for dir in xvg_dirs:
            os.makedirs(dir, exist_ok=True)
        
        run_isomerization_analysis(args, scripts_dir, out_dir)

    elif args.analysis == "process":
        run_trajectory_processing(args, scripts_dir, out_dir, tpr_path, xtc_path, gro_path)
    elif args.analysis == "msd":
        xvg_dir = os.path.join(out_dir, "xvg_files")
        os.makedirs(xvg_dir, exist_ok=True)
        run_msd_analysis(args, scripts_dir, out_dir, tpr_path, xtc_path)
    elif args.analysis == "diffusion":
        run_diffusion_analysis(args, scripts_dir, out_dir, xvg_path, outputs_dir)



def run_trajectory_processing(args, scripts_dir, out_dir, tpr_path, xtc_path, gro_path):
    script_path = os.path.join(scripts_dir, "process_trajectory.py")
    cmd = [sys.executable, script_path]
    
    if not tpr_path:
        print("tpr required")
        sys.exit(1)
    cmd.extend(["--tpr", tpr_path])
    if not (gro_path or xtc_path):
        print("xtc or gro required")
        sys.exit(1)
    if gro_path:
        cmd.extend(["--gro", gro_path])
    if xtc_path:
        cmd.extend(["--xtc", xtc_path])
    if args.begin_time is not None:
        cmd.extend(["--begin-time", str(args.begin_time)])
    if args.skip is not None:
        cmd.extend(["--skip", str(args.skip)])
    if args.lipid_selection:
        cmd.extend(["--lipid-selection", args.lipid_selection])
    elif args.resname:
        cmd.extend(["--resname", args.resname])
    cmd.extend(["--out_dir", out_dir])
    
    print(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def run_dihedral_analysis(args, scripts_dir, out_dir, xvg_dir, idx_dir, xtc_path):
    script_path = os.path.join(scripts_dir, "dihedral_distributions.py")
    cmd = [sys.executable, script_path]
    
    if not xtc_path:
        print("xtc required")
        sys.exit(1)
    cmd.extend(["--xtc", xtc_path])
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
    subprocess.run(cmd, check=True)

def run_angle_analysis(args, scripts_dir, out_dir, xvg_dir, idx_dir, xtc_path):
    script_path = os.path.join(scripts_dir, "angle_distributions.py")
    cmd = [sys.executable, script_path]
    if not xtc_path:
        print("xtc required")
        sys.exit(1)
    cmd.extend(["--xtc", xtc_path])
    cmd.extend(["--out_dir", xvg_dir])
    cmd.extend(["--idx_dir", idx_dir])
    if args.nlip:
        cmd.extend(["--nlip", str(args.nlip)])
    if args.offset:
        cmd.extend(["--offset", str(args.offset)])
    if args.angle_set:
        cmd.extend(["--angle_set", args.angle_set])
    if args.xvg_suffix:
        cmd.extend(["--xvg-suffix", args.xvg_suffix])
    
    print(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def run_apl_analysis(args, scripts_dir, out_dir, gro_path, xtc_path):
    script_path = os.path.join(scripts_dir, "APLfromTRJ.py")
    cmd = [sys.executable, script_path]
    if not gro_path:
        print("gro required")
        sys.exit(1)
    cmd.extend(["--gro", gro_path])
    if not xtc_path:
        print("xtc required")
        sys.exit(1)
    cmd.extend(["--xtc", xtc_path])
    if not (args.resname or args.lipid21_resnames):
        print("Resname required")
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
    cmd.extend(["--out_dir", out_dir])
    
    print(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def run_rdf_analysis(args, scripts_dir, out_dir, xtc_path, tpr_path):
    script_path = os.path.join(scripts_dir, "RDF_multiprocess_PO-NC_IDexclusion.py")
    cmd = [sys.executable, script_path]

    if not xtc_path:
        print("xtc required")
        sys.exit(1)
    cmd.extend(["--xtc", xtc_path])
    if not tpr_path:
        print("tpr required")
        sys.exit(1)
    cmd.extend(["--tpr", tpr_path])
    cmd.extend(["--out_dir", out_dir])
    index_dir = os.path.join(out_dir, "index_files")
    xvg_dir = os.path.join(out_dir, "xvg_files")
    cmd.extend(["--index-dir", index_dir])
    cmd.extend(["--xvg-dir", xvg_dir])
    
    print(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def run_isomerization_analysis(args, scripts_dir, out_dir):
    script_path = os.path.join(scripts_dir, "isomerization_code_PointGraph.py")
    cmd = [sys.executable, script_path]
    
    data_dir = out_dir
    os.makedirs(data_dir, exist_ok=True)
    
    # this is temporary data. will add lipid21 and charmm36
    xvg_dirs = [
        os.path.join(data_dir, "xvg_files_OFF_HMR"),
        os.path.join(data_dir, "xvg_files_macrog"), 
        os.path.join(data_dir, "xvg_files_slipids")
    ]
    for dir in xvg_dirs:
        os.makedirs(dir, exist_ok=True)
    output_file = os.path.join(out_dir, "plot_outputs", "sn1_isomerization_comparison.png")
    cmd.extend(["--xvg_dirs"] + xvg_dirs)
    cmd.extend(["--out_dir", out_dir])
    cmd.extend(["--output", output_file])
    if args.title:
        cmd.extend(["--title", args.title])
    
    print(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def run_msd_analysis(args, scripts_dir, out_dir, tpr_path, xtc_path):
    script_path = os.path.join(scripts_dir, "msd_plot.py")
    cmd = [sys.executable, script_path]
    
    if not tpr_path:
        print("No tpr")
        sys.exit(1)
    cmd.extend(["--tpr", tpr_path])
    if xtc_path:
        cmd.extend(["--xtc", xtc_path])
    if args.prefix:
        cmd.extend(["--prefix", args.prefix])
    cmd.extend(["--phosphate-atom", args.phosphate_atom])
    cmd.extend(["--out_dir", out_dir])
    
    print(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

# added this after a hiatus and did not revisit structure
def run_diffusion_analysis(args, scripts_dir, out_dir, xvg_path, outputs_dir):
    script_path = os.path.join(scripts_dir, "lateral_diffusion_quant.py")
    cmd = [sys.executable, script_path]
    
    #manual default for xvg file made by previous msd function
    if xvg_path:
        cmd.extend(["--xvg", xvg_path])
    else:
        default_msd_xvg_dir = os.path.join(outputs_dir, "msd_output", "xvg_files") 
        if os.path.exists(default_msd_xvg_dir):
            xvg_files = [f for f in os.listdir(default_msd_xvg_dir) if f.endswith('.xvg')]
            if xvg_files:
                default_xvg_path = os.path.join(default_msd_xvg_dir, xvg_files[0])
                print(f"Using default xvg from msd_plot.py: {default_xvg_path}")
                cmd.extend(["--xvg", default_xvg_path])
            else:
                print(f"No defualt xvgs: {default_msd_xvg_dir}")
        else:
            print(f"No default directory: {default_msd_xvg_dir}")

    if args.window_size:
        cmd.extend(["--window-size", str(args.window_size)])
    if args.window_step:
        cmd.extend(["--window-step", str(args.window_step)])
    if args.min_r2:
        cmd.extend(["--min-r2", str(args.min_r2)])
    if args.title:
        cmd.extend(["--title", args.title])
    cmd.extend(["--out_dir", out_dir])
    
    print(f"Running command: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

if __name__ == "__main__":
    main()