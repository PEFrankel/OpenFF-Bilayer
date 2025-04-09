import argparse
import subprocess
import sys

def run_command(script, args_list, stats):
    print(f"\n===== {script} =====")
    print(f"Running: {' '.join(args_list)}")
    try:
        subprocess.run(args_list, check=True)
        stats["passed"] += 1
    except subprocess.CalledProcessError as e:
        stats["failed"] += 1

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--xtc", required=True, help="Path to XTC file")
    parser.add_argument("--tpr", required=True, help="Path to TPR file")
    parser.add_argument("--gro", required=True, help="Path to GRO file")
    parser.add_argument("--resname", required=True, help="Lipid residue name")
    parser.add_argument("--nlip", type=int, required=True, help="Number of lipids")
    parser.add_argument("--cli-path", default="./OFFCLI.py")
    args = parser.parse_args()

    py = sys.executable
    cli = args.cli_path

    stats = {"passed": 0, "failed": 0}

    tests = [
    #processing is modified to begin at 100 ns and skip 10 ps (much faster and the extra data does not alter results significantly)
    ("Testing process",[py, cli, "--analysis", "process", "--xtc", args.xtc, "--tpr", args.tpr, "--gro", args.gro, "--resname", args.resname,"--begin-time", "100000", "--skip", "10"]),
    ("Testing dihedral",[py, cli, "--analysis", "dihedral", "--xtc", args.xtc, "--torsion_set", "sn1_isomerization", "--nlip", str(args.nlip)]),
    ("Testing angle", [py, cli, "--analysis", "angle", "--xtc", args.xtc, "--angle_set", "phosphate", "--nlip", str(args.nlip)]),
    ("Testing apl", [py, cli, "--analysis", "apl", "--gro", args.gro, "--xtc", args.xtc, "--resname", args.resname]),
    ("Testing rdf", [py, cli, "--analysis", "rdf", "--xtc", args.xtc, "--tpr", args.tpr]),
    ("Testing msd", [py, cli, "--analysis", "msd", "--xtc", args.xtc, "--tpr", args.tpr]),
    # window size and step are equivalent to 10 picoseconds i.e. window size of 10000 and step of 200 covers 100 ns of 2 ns increments
    ("Testing diffusion", [py, cli, "--analysis", "diffusion", "--window-size", "1000", "--window-step", "20"]), # this is -skip 10 specific (ratio of 50 for good 100ns)
    # isomerization code doesn't read into the dihedral output rn. nRFT
    # ("Testing isomerization", [py, cli, "--analysis", "isomerization"]),
    ]

    for desc, cmd in tests:
        run_command(desc, cmd, stats)

    print("\n \_DONE_/")
    print(f"Passed: {stats['passed']}")
    print(f"Failed: {stats['failed']}")

if __name__ == "__main__":
    main()
