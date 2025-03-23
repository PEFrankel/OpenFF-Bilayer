import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
import argparse
import os
import sys

def calculate_area_per_lipid(universe, lipid_resnames, standard_resname):
    if standard_resname:
        # For standard POPC or specified standard residue name
        lipid_group = universe.select_atoms(f"resname {standard_resname}")
        n_lipids = len(lipid_group.residues)
    else:
        # For files where POPC is split into multiple residues (Lipid17/21)
        lipid_group = universe.select_atoms(f"resname {' '.join(lipid_resnames)}")
        n_lipids = len(lipid_group.residues) // len(lipid_resnames) if lipid_group and len(lipid_resnames) > 0 else 0

    print(f"Number of lipids detected: {n_lipids}")
    apl = {}

    for ts in universe.trajectory:
        dimensions = ts.dimensions
        apl_frame = (dimensions[0] * dimensions[1] * 2) / n_lipids if n_lipids > 0 else np.nan  # APL calc with box area
        apl[ts.time] = apl_frame

    return apl, n_lipids

def plot_area_per_lipid(trajectory_file, structure_file, lipid21_resnames, standard_resname, title, color, max_time=None, output_prefix=None):
    plt.figure(figsize=(10, 6))

    u = mda.Universe(structure_file, trajectory_file)
    
    apl_data, n_lipids = calculate_area_per_lipid(u, lipid21_resnames, standard_resname)
    times = np.array(list(apl_data.keys()))
    area_per_lipid = np.array(list(apl_data.values()))

    # Apply time filter if specified
    if max_time:
        filtered_indices = np.where(times <= max_time)[0] # Within specified filter
        filtered_times = times[filtered_indices]
        filtered_area_per_lipid = area_per_lipid[filtered_indices]
    else:
        filtered_times = times
        filtered_area_per_lipid = area_per_lipid

    # Calculate statistics
    average_apl = np.nanmean(filtered_area_per_lipid)
    std_apl = np.nanstd(filtered_area_per_lipid)
    n_frames = len(filtered_times)  # Standard error over total ensemble frames
    sem_apl = std_apl / np.sqrt(n_frames)
    
    print(f"Average APL: {average_apl:.2f} Å²")
    print(f"Standard Deviation: {std_apl:.2f} Å²")
    print(f"Standard Error of Mean: {sem_apl:.4f} Å²")
    print(f"Number of frames: {n_frames}")

    # Plot
    filtered_times_ns = filtered_times / 1000
    plt.plot(filtered_times_ns, filtered_area_per_lipid, alpha=0.5, 
             label=f"{title}: {average_apl:.2f} ± {sem_apl:.2f} Å²", 
             color=color)
    
    # Calculate and display trend line (disabled)
    coefficients = np.polyfit(filtered_times, filtered_area_per_lipid, 1)
    x_fit = np.linspace(min(filtered_times), max(filtered_times), 1000)
    y_fit = np.polyval(coefficients, x_fit)
    x_fit_ns = x_fit / 1000  # Convert to nanoseconds
    # plt.plot(x_fit_ns, y_fit, '--', color=color, alpha=0.7)

    print(f"Linear fit - Slope: {coefficients[0]:.6f}, Intercept: {coefficients[1]:.6f}")

    plt.xlabel('Time (ns)')
    plt.ylabel('Area per lipid (Å²)')
    plt.title(title)
    # plt.xlim[(0, 500)]
    plt.ylim([45, 75])
    plt.legend()
    
    # Set output filename with prefix if provided
    if output_prefix:
        output_filename = f"{output_prefix}"
    else:
        output_filename = f"{title.replace(' ', '_')}_APL"
    
    # Save figure
    plt.savefig(f"{output_filename}.png", dpi=300)
    
    # Data saves to text file if interested
    with open(f"{output_filename}.txt", 'w') as f:
        f.write(f"# {title} - Area per lipid analysis\n")
        f.write(f"# Average APL: {average_apl:.4f} Å²\n")
        f.write(f"# Standard Deviation: {std_apl:.4f} Å²\n")
        f.write(f"# Standard Error of Mean: {sem_apl:.4f} Å² (calculated using {n_frames} frames)\n")
        f.write(f"# Linear fit - Slope: {coefficients[0]:.6f}, Intercept: {coefficients[1]:.6f}\n")
        f.write("# Time(ps) Area_per_lipid(Å²)\n")
        np.savetxt(f, np.column_stack((times, area_per_lipid)))
    
    return output_filename

def main():
    # Color options with corresponding force fields
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

    parser = argparse.ArgumentParser(description="Calculate and analyze area per lipid from trajectory")
    parser.add_argument("--gro", required=True, help="Path to GRO structure file")
    parser.add_argument("--xtc", required=True, help="Path to XTC trajectory file")
    parser.add_argument("--title", default="Area Per Lipid Analysis", 
                        help="Title for the plot and output files")
    parser.add_argument("--max-time", type=int, default=1000000, 
                        help="Maximum time (ps) to include in analysis (default: 1000000)")
    parser.add_argument("--color", default="openff", choices=list(color_choices.keys()),
                        help="Color to use for plotting. Options include force field identifiers (openff, slipids, lipid21, macrog, charmm36, core, aux) or basic colors.")
    parser.add_argument("--output", default=None, help="Output file prefix (without extension)")
    
    # Mutually exclusive group for residue naming options
    residue_group = parser.add_mutually_exclusive_group(required=True)
    residue_group.add_argument("--standard-resname", help="Standard residue name (e.g., POPC)")
    residue_group.add_argument("--lipid21-resnames", nargs='+', help="Residue names for split lipids (e.g., PA PC OL for Lipid17/21)")
    
    args = parser.parse_args()

    # Validate input files
    if not os.path.exists(args.gro):
        print(f"Error: GRO file '{args.gro}' does not exist")
        sys.exit(1)
    if not os.path.exists(args.xtc):
        print(f"Error: XTC file '{args.xtc}' does not exist")
        sys.exit(1)
    
    # Get color from choices dict
    color = color_choices.get(args.color.lower(), "#0333b0")  # Default to blue if not found
    
    # Create output directory if needed
    if args.output:
        output_dir = os.path.dirname(args.output)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
    
    # Run analysis
    output_name = plot_area_per_lipid(
        trajectory_file=args.xtc,
        structure_file=args.gro,
        lipid21_resnames=args.lipid21_resnames if args.lipid21_resnames else [],
        standard_resname=args.standard_resname,
        title=args.title,
        color=color,
        max_time=args.max_time,
        output_prefix=args.output
    )
    
    print(f"Analysis complete. Results saved to {output_name}.png and {output_name}.txt")

if __name__ == "__main__":
    main()