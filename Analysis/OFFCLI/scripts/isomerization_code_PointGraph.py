import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import MDAnalysis as mda
import argparse

def main():
    parser = argparse.ArgumentParser(description="Calculate sn-1 gauche isomerization from dihedral angle distributions")
    parser.add_argument("--xvg_folder", help="Folder containing XVG files (default: uses predefined list)", default=None)
    parser.add_argument("--xvg_folders", nargs='+', help="Multiple folders containing XVG files", default=None)
    parser.add_argument("--output", help="Output file path for the plot (default: sn1_isomerization_comparison.png)", 
                        default="sn1_isomerization_comparison.png")
    parser.add_argument("--title", help="Plot title (default: sn-1 Gauche Isomerization)", 
                        default="sn-1 Gauche Isomerization")
    parser.add_argument("--positive_range", nargs=2, type=float, default=(30, 90),
                        help="Positive angle range (default: 30 90)")
    parser.add_argument("--negative_range", nargs=2, type=float, default=(-90, -30),
                        help="Negative angle range (default: -90 -30)")
    parser.add_argument("--color_scheme", nargs='+', default=['#0333b0', 'green', '#ee7f17'],
                        help="Color scheme for plots (default: '#0333b0' 'green' '#ee7f17')")
    parser.add_argument("--labels", nargs='+', default=None,
                        help="Labels for each folder (default: folder names)")
    
    args = parser.parse_args()
    
    # Set up folders with xvg files
    if args.xvg_folders is not None:
        xvg_folders = args.xvg_folders
    elif args.xvg_folder is not None:
        xvg_folders = [args.xvg_folder]
    else:
        # Use default folders if none provided
        xvg_folders = ['xvg_files_OFF_HMR', 'xvg_files_macrog', 'xvg_files_slipids']
    
    # Validate that folders exist
    for folder in xvg_folders:
        if not os.path.exists(folder):
            print(f"Warning: XVG folder '{folder}' does not exist. Skipping.")
            xvg_folders.remove(folder)
    
    if len(xvg_folders) == 0:
        print("Error: No valid XVG folders found. Exiting.")
        return 1
    
    # Set up labels
    labels = args.labels
    if labels is None:
        # Use folder names as labels by default
        labels = [os.path.basename(folder).replace('xvg_files_', '') for folder in xvg_folders]
        # Ensure we have enough labels
        if len(labels) < len(xvg_folders):
            labels.extend([f"Dataset {i+1}" for i in range(len(labels), len(xvg_folders))])
    
    # Ensure we have enough colors
    colors = args.color_scheme
    if len(colors) < len(xvg_folders):
        # Add more colors if needed
        default_colors = ['blue', 'green', 'red', 'purple', 'orange', 'brown', 'pink', 'gray', 'olive', 'cyan']
        colors.extend([default_colors[i % len(default_colors)] for i in range(len(colors), len(xvg_folders))])
    
    # Set angle ranges
    positive_range = tuple(args.positive_range)
    negative_range = tuple(args.negative_range)
    
    # Dictionary to store percentages for each folder
    percentages_dict = {folder: [] for folder in xvg_folders}
    
    # Calculate probability percentages within ranges for each XVG folder
    for folder in xvg_folders:
        xvg_files = [os.path.join(folder, f) for f in sorted(os.listdir(folder)) if f.endswith('.xvg')]
        print(f"Processing {len(xvg_files)} XVG files in {folder}")
        
        for xvg_file in xvg_files:
            try:
                combined_percent = calculate_percentage(xvg_file, positive_range, negative_range)
                percentages_dict[folder].append(combined_percent)
            except Exception as e:
                print(f"Error processing {xvg_file}: {e}")
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10, 6))
    
    for folder, color, label in zip(xvg_folders, colors, labels):
        percentages = percentages_dict[folder]
        if not percentages:
            print(f"Warning: No valid data for {folder}, skipping in plot")
            continue
            
        x_labels = np.arange(1, len(percentages) + 1)
        mean_percentage = np.mean(percentages)
        se_percentage = np.std(percentages) / np.sqrt(len(percentages)) if len(percentages) > 1 else 0
        
        ax.errorbar(x_labels, percentages, yerr=se_percentage, fmt='o', color=color, 
                   ecolor=color, capsize=5, 
                   label=f'{label} (Mean: {mean_percentage:.2f}%, SE: {se_percentage:.2f}%)')
    
    ax.set_xlabel('Bond #')
    ax.set_ylabel('% Isomerization')
    ax.set_title(args.title)
    ax.legend()
    
    # Set x-ticks based on the maximum number of bonds across all datasets
    max_bonds = max([len(percentages_dict[folder]) for folder in xvg_folders])
    ax.set_xticks(np.arange(1, max_bonds + 1))
    
    fig.tight_layout()
    
    # Save the plot
    plt.savefig(args.output, dpi=300)
    print(f"Plot saved to {args.output}")
    
    # Output the results into XVG file for each folder
    for folder in xvg_folders:
        percentages = percentages_dict[folder]
        if not percentages:
            continue
            
        output_xvg = f'sn1_isomerization_{os.path.basename(folder)}.xvg'
        with open(output_xvg, 'w') as file:
            file.write('# Probability Percentages within Specified Angle Ranges\n')
            file.write(f'# Bond#  Combined({positive_range[0]} to {positive_range[1]} and {negative_range[0]} to {negative_range[1]})\n')
            for i, combined_percent in enumerate(percentages, start=1):
                file.write(f"{i}  {combined_percent:.2f}\n")
        print(f"Results saved to {output_xvg}")
    
    # Optionally display the plot
    plt.show()
    
    return 0

def calculate_percentage(xvg_file, positive_range, negative_range):
    """Calculate probability percentage within specified angle ranges from an XVG file"""
    try:
        aux = mda.auxiliary.XVG.XVGReader(xvg_file)
        data = []
        for step in aux:
            step_data = step.data
            if step_data.ndim == 1:
                step_data = step_data.reshape(-1, 2)
            data.extend(step_data)
        
        df = pd.DataFrame(data, columns=['A', 'B'])
        
        # Sum of probabilities within ranges
        positive_prob = df[(df['A'] >= positive_range[0]) & (df['A'] <= positive_range[1])]['B'].sum()
        negative_prob = df[(df['A'] >= negative_range[0]) & (df['A'] <= negative_range[1])]['B'].sum()
        
        # Calculate total probability to normalize
        total_prob = df['B'].sum()
        
        # Sum
        combined_prob = positive_prob + negative_prob
        
        # Convert sum to percentage
        combined_percentage = (combined_prob / total_prob) * 100 if total_prob != 0 else 0
        
        return combined_percentage
    except Exception as e:
        print(f"Error processing file {xvg_file}: {e}")
        raise

if __name__ == "__main__":
    exit_code = main()
    exit(exit_code)