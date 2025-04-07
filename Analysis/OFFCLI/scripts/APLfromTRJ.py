import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
import argparse
import os
import sys
import json

def calculate_area_per_lipid(universe, lipid_resnames, standard_resname):
    if standard_resname:
        lipid_group = universe.select_atoms(f"resname {standard_resname}")
        n_lipids = len(lipid_group.residues)
    else:
        # for lipid21 (split residues)
        lipid_group = universe.select_atoms(f"resname {' '.join(lipid_resnames)}")
        n_lipids = len(lipid_group.residues) // len(lipid_resnames) if lipid_group and len(lipid_resnames) > 0 else 0

    # print(f"Number of lipids: {n_lipids}")
    apl = {}

    for ts in universe.trajectory:
        dimensions = ts.dimensions
        apl_frame = (dimensions[0] * dimensions[1] * 2) / n_lipids if n_lipids > 0 else np.nan  # APL calc with box area
        apl[ts.time] = apl_frame

    return apl, n_lipids

def plot_area_per_lipid(trajectory_file, structure_file, lipid21_resnames, standard_resname, title, color, out_dir, max_time=None, output_prefix=None):
    os.makedirs(out_dir, exist_ok=True)
    plot_outputs_dir = os.path.join(out_dir, "plot_outputs")
    os.makedirs(plot_outputs_dir, exist_ok=True)

    plt.figure(figsize=(12,8))

    u = mda.Universe(structure_file, trajectory_file)
    
    apl_data, n_lipids = calculate_area_per_lipid(u, lipid21_resnames, standard_resname)
    times = np.array(list(apl_data.keys()))
    area_per_lipid = np.array(list(apl_data.values()))

    #filter, probably unused
    if max_time:
        filtered_indices = np.where(times <= max_time)[0]
        filtered_times = times[filtered_indices]
        filtered_area_per_lipid = area_per_lipid[filtered_indices]
    else:
        filtered_times = times
        filtered_area_per_lipid = area_per_lipid

    average_apl = np.nanmean(filtered_area_per_lipid)
    std_apl = np.nanstd(filtered_area_per_lipid)
    n_frames = len(filtered_times)               # SE over frames in sims
    sem_apl = std_apl / np.sqrt(n_frames)
    
    # print(f"Avg APL: {average_apl:.2f} Å²")
    # print(f"STD: {std_apl:.2f} Å²")
    # print(f"SEM: {sem_apl:.4f} Å²")
    # print(f"Number of frames: {n_frames}")

    filtered_times_ns = filtered_times / 1000
    plt.plot(filtered_times_ns, filtered_area_per_lipid, alpha=0.5, 
             label=f"{title}: {average_apl:.2f} ± {sem_apl:.2f} Å²", 
             color=color)
    
    # Calculate and display trend line (DISABLED)
    coefficients = np.polyfit(filtered_times, filtered_area_per_lipid, 1)
    x_fit = np.linspace(min(filtered_times), max(filtered_times), 1000)
    y_fit = np.polyval(coefficients, x_fit)
    x_fit_ns = x_fit / 1000  # ns
    # plt.plot(x_fit_ns, y_fit, '--', color=color, alpha=0.7)

    print(f"Linear fit - Slope: {coefficients[0]:.6f}, Intercept: {coefficients[1]:.6f}")

    plt.xlabel('Time (ns)')
    plt.ylabel('Area per lipid (Å²)')
    plt.title(title)
    # plt.xlim[(0, 500)]
    plt.ylim([45, 75]) # adjust this for continuity with other ff plots
    plt.legend()
    
    if output_prefix:
        output_filename = os.path.join(plot_outputs_dir, f"{output_prefix}")
    else:
        output_filename = os.path.join(plot_outputs_dir, f"{title.replace(' ', '_')}_APL")
    
    plt.tight_layout()
    plt.savefig(f"{output_filename}.png", dpi=300)
    
    with open(f"{output_filename}.txt", 'w') as f: #data wrote to txt file. for any post process
        f.write(f"# {title} - Area per lipid analysis\n")
        f.write(f"# Average APL: {average_apl:.4f} Å²\n")
        f.write(f"# STD: {std_apl:.4f} Å²\n")
        f.write(f"# SEM: {sem_apl:.4f} Å² (calculated using {n_frames} frames)\n")
        # f.write(f"# Linear fit - Slope: {coefficients[0]:.6f}, Intercept: {coefficients[1]:.6f}\n")
        f.write("# Time(ps) Area_per_lipid(Å²)\n")
        np.savetxt(f, np.column_stack((times, area_per_lipid)))
    
        json_output = {
        "title": title,
        "average_apl": float(average_apl),
        "standard_deviation": float(std_apl),
        "standard_error_of_mean": float(sem_apl),
        "number_of_frames": int(n_frames),
        "linear_fit": {
            "slope": float(coefficients[0]),
            "intercept": float(coefficients[1])
        }
    }

    if output_prefix:
        json_filename = os.path.join(plot_outputs_dir, f"{output_prefix}_stats.json")
    else:
        json_filename = os.path.join(plot_outputs_dir, f"{title.replace(' ', '_')}_APL_stats.json")
    
    with open(json_filename, 'w') as f:
        json.dump(json_output, f, indent=4)

    return os.path.basename(output_filename)

def main():
    color_choices = {
        "openff": "#0333b0",    # OpenFF
        "slipids": "#ee7f17",   # Slipids
        "lipid21": "#ee7f17",   # Lipid21
        "macrog": "green",      # MacRog
        "charmm36": "green",    # CHARMM36
        "core": "black",        # core
        "aux": "purple",        # aux
        "blue": "blue",
        "red": "red",
        "orange": "orange",
        "green": "green",
        "black": "black",
        "purple": "purple",
        "grey": "grey"
    }

    parser = argparse.ArgumentParser(description="Calculate and plot area per lipid")
    parser.add_argument("--gro", required=True, help="Path to GRO file")
    parser.add_argument("--xtc", required=True, help="Path to XTC file")
    parser.add_argument("--out_dir", required=True, help="Output directory")
    parser.add_argument("--title", default="Area Per Lipid Analysis", help="Title for the plot and output files")
    parser.add_argument("--max-time", type=int, default=1000000, help="Maximum time (ps) to include")
    parser.add_argument("--color", default="openff", choices=list(color_choices.keys()),
    help="Options for FFs (openff, slipids, lipid21, macrog, charmm36, core, aux) or basic colors.")
    
    residue_group = parser.add_mutually_exclusive_group(required=True)
    residue_group.add_argument("--standard-resname", help="Standard residue name (POPC)")
    residue_group.add_argument("--lipid21-resnames", nargs='+', help="Residue names for split lipids (e.g. PA PC OL for Lipid17/21)")
    args = parser.parse_args()
    
    color = color_choices.get(args.color.lower(), "#0333b0")  # OFF default
    
    output_name = plot_area_per_lipid(
        trajectory_file=args.xtc,
        structure_file=args.gro,
        lipid21_resnames=args.lipid21_resnames if args.lipid21_resnames else [],
        standard_resname=args.standard_resname,
        title=args.title,
        color=color,
        out_dir=args.out_dir,
        max_time=args.max_time
    )
    
    # print(f"Done, saved to {output_name}")

if __name__ == "__main__":
    main()