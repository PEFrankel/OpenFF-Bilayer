import argparse
import os
import sys
import tempfile
import subprocess


def process_trajectory(tpr_file, gro_file=None, xtc_file=None, output_dir="./processed",begin_time=None, skip=None, lipid_selection=None):
    # if not os.path.exists(tpr_file):
    #     raise FileNotFoundError(f"no tpr: {tpr_file}")
    # if gro_file is None and xtc_file is None:
    #     raise ValueError("no gro or xtc")
    # if gro_file and not os.path.exists(gro_file):
    #     raise FileNotFoundError(f"no gro: {gro_file}")
    # if xtc_file and not os.path.exists(xtc_file):
    #     raise FileNotFoundError(f"no xtc: {xtc_file}")
    ##^^debugs^^##

    os.makedirs(output_dir, exist_ok=True)
    data_dir = os.path.join(output_dir, "data")
    os.makedirs(data_dir, exist_ok=True)
    tpr_file = os.path.abspath(tpr_file)
    if gro_file:
        gro_file = os.path.abspath(gro_file)
    if xtc_file:
        xtc_file = os.path.abspath(xtc_file)
    
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as system_sel:
        system_sel.write("System\n")
        system_sel_path = system_sel.name

    if not lipid_selection:
        lipid_selection = "Lipid"
    
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as lipid_sel:
        lipid_sel.write(f"{lipid_selection}\n")
        lipid_sel.write("System\n")
        lipid_sel_path = lipid_sel.name
    
    output_files = {"gro": None, "xtc": None}
    
    if gro_file:
        nojump_gro = os.path.join(data_dir, "nojump.gro")
        # print(f"GMX COMMAND: gmx trjconv -s {tpr_file} -f {gro_file} -pbc nojump -o {nojump_gro} <<<<<<<<<<<")
        
        subprocess.run(
            f"gmx trjconv -s {tpr_file} -f {gro_file} -pbc nojump -o {nojump_gro}",
            shell=True, 
            stdin=open(system_sel_path, 'r'),
            check=True
        )
        
        centered_gro = os.path.join(data_dir, "centered.gro")
        subprocess.run(
            f"gmx trjconv -s {tpr_file} -f {nojump_gro} -pbc mol -center -o {centered_gro}",
            shell=True,
            stdin=open(lipid_sel_path, 'r'),
            check=True
        )
        
        analysis_gro = os.path.join(data_dir, "analysis.gro")
        subprocess.run(
            f"gmx trjconv -s {tpr_file} -f {centered_gro} -pbc whole -o {analysis_gro}",
            shell=True,
            stdin=open(system_sel_path, 'r'),
            check=True
        )
        
        output_files["gro"] = analysis_gro
        # print(f"gro output dir: {analysis_gro}")
        
        try:
            os.remove(nojump_gro)
            os.remove(centered_gro)
        except OSError as e:
            print(f"temps: {e}")
    
    # now TRJ
    if xtc_file:
        time_param = f"-b {begin_time}" if begin_time is not None else ""
        skip_param = f"-skip {skip}" if skip is not None else ""
        nojump_xtc = os.path.join(data_dir, "nojump.xtc")

        # print(f"GMX COMMAND: gmx trjconv -s {tpr_file} -f {xtc_file} {time_param} {skip_param} -pbc nojump -o {nojump_xtc} <<<<<<<<<<<")
        subprocess.run(
            f"gmx trjconv -s {tpr_file} -f {xtc_file} {time_param} {skip_param} -pbc nojump -o {nojump_xtc}",
            shell=True,
            stdin=open(system_sel_path, 'r'),
            check=True
        )
        
        centered_xtc = os.path.join(data_dir, "centered.xtc")
        subprocess.run(
            f"gmx trjconv -s {tpr_file} -f {nojump_xtc} -pbc mol -center -o {centered_xtc}",
            shell=True,
            stdin=open(lipid_sel_path, 'r'),
            check=True
        )
        
        analysis_xtc = os.path.join(data_dir, "analysis.xtc")
        subprocess.run(
            f"gmx trjconv -s {tpr_file} -f {centered_xtc} -pbc whole -o {analysis_xtc}",
            shell=True,
            stdin=open(system_sel_path, 'r'),
            check=True
        )
        
        output_files["xtc"] = analysis_xtc
        # print(f"xtc output: {analysis_xtc}")
        
        try:
            os.remove(nojump_xtc)
            os.remove(centered_xtc)
        except OSError as e:
            print(f"temps: {e}")
    
    os.unlink(system_sel_path)
    os.unlink(lipid_sel_path)
    
    return output_files


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.dirname(script_dir)
    
    parser = argparse.ArgumentParser(description="Process GROMACS trajectory files for analysis")
    parser.add_argument("--tpr", required=True, help="Path to TPR file")
    parser.add_argument("--gro", help="Path to GRO file")
    parser.add_argument("--xtc", help="Path to XTC file")
    parser.add_argument("--out_dir", default=os.path.join(base_dir, "outputs", "process_output"), help="dir = outputs/process_output")
    parser.add_argument("--begin-time", type=int, help="Start time (ps) for trj processing")
    parser.add_argument("--skip", type=int, help="Frame skip for trj processing")
    parser.add_argument("--resname", help="Standard residue name for lipid selection (POPC)")
    parser.add_argument("--lipid-selection", help="Alt selection for lipids (default: 'Lipid' or resname)")
    args = parser.parse_args()
    
    # if not args.gro and not args.xtc:
    #     parser.error("no gro and/or xtc")
    
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
            output_dir=args.out_dir,
            begin_time=args.begin_time,
            skip=args.skip,
            lipid_selection=lipid_selection
        )
        #Result check
        # print("\nTrajectory processing completed successfully.")
        # print(f"gro: {output_files['gro'] if output_files['gro'] else 'fail'}")
        # print(f"xtc: {output_files['xtc'] if output_files['xtc'] else 'fail'}")
        
    except Exception as e:
        print(f"error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()