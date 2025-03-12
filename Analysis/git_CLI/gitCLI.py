import argparse
import os
import subprocess
import sys

def main():
    parser = argparse.ArgumentParser(description="Lipid bilayer analysis tool")
    parser.add_argument("--analysis", choices=["dihedral", "angle", "apl", "rdf", "isomerization"], 
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

    args = parser.parse_args()

    # Run the appropriate script based on the analysis type
    if args.analysis == "dihedral":
        run_dihedral_analysis(args)
    elif args.analysis == "angle":
        run_angle_analysis(args)
    elif args.analysis == "apl":
        run_apl_analysis(args)
    elif args.analysis == "rdf":
        run_rdf_analysis(args)
    elif args.analysis == "isomerization":
        run_isomerization_analysis(args)

def run_dihedral_analysis(args):
    cmd = [sys.executable, "dihedral_distributions.py"]
    
    if not args.xtc:
        print("Error: --xtc is required for dihedral analysis")
        sys.exit(1)
    
    cmd.extend(["--xtc", args.xtc])
    
    if args.nlip:
        cmd.extend(["--nlip", str(args.nlip)])
    if args.offset:
        cmd.extend(["--offset", str(args.offset)])
    if args.torsion_set:
        cmd.extend(["--torsion_set", args.torsion_set])
    
    subprocess.run(cmd)

def run_angle_analysis(args):
    cmd = [sys.executable, "angle_distributions.py"]
    
    # environment variables
    env = os.environ.copy()
    if args.xtc:
        env["TRAJECTORY_FILE"] = args.xtc
    if args.nlip:
        env["NLIP"] = str(args.nlip)
    if args.offset:
        env["OFFSET"] = str(args.offset)
    
    subprocess.run(cmd, env=env)

def run_apl_analysis(args):
    cmd = [sys.executable, "APLfromTRJ.py"]
    
    if not args.gro:
        print("Error: --gro is required for APL analysis")
        sys.exit(1)
    cmd.extend(["--gro", args.gro])
    if not args.xtc:
        print("Error: --xtc is required for APL analysis")
        sys.exit(1)
    cmd.extend(["--xtc", args.xtc])
    
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
    
    subprocess.run(cmd)

def run_rdf_analysis(args):
    cmd = [sys.executable, "RDF_multiprocess_PO-NC_IDexclusion.py"]
    
    # environment variables
    env = os.environ.copy()
    if args.xtc:
        env["XTC_FILE"] = args.xtc
    if args.tpr:
        env["TPR_FILE"] = args.tpr
    
    subprocess.run(cmd, env=env)

def run_isomerization_analysis(args):
    cmd = [sys.executable, "isomerization_code_PointGraph.py"]
    
    # environment variables
    env = os.environ.copy()
    
    subprocess.run(cmd, env=env)

if __name__ == "__main__":
    main()