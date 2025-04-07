import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import MDAnalysis as mda
import argparse

def main():
    parser = argparse.ArgumentParser(description="Calculate gauche isomerization from sn-1 torsion distributions")
    parser.add_argument("--xvg_dir", help="Path to XVG files", default=None)#debug, used >>>>>>>>>remove these duplicates
    parser.add_argument("--xvg_dirs", nargs='+', help="Path to XVG files", default=None) #same as above
    parser.add_argument("--output", help="Output for the plot (default: plot_outputs/sn1_isomerization_comparison.png)", default=None)#debug, still usable
    parser.add_argument("--out_dir", help="Output directory", required=True) #again, effectively same as above. 
    parser.add_argument("--title", help="Plot title (default: sn-1 Gauche Isomerization)", default="sn-1 Gauche Isomerization")
    parser.add_argument("--positive_range", nargs=2, type=float, default=(30, 90),help="Default: 30 90")
    parser.add_argument("--negative_range", nargs=2, type=float, default=(-90, -30), help="Default: -90 -30")
    parser.add_argument("--color_scheme", nargs='+', default=['#0333b0', 'green', '#ee7f17'], help="Color scheme for plot (default: '#0333b0' 'green' '#ee7f17')")
    parser.add_argument("--labels", nargs='+', default=None, help="Default: dir names")
    args = parser.parse_args()
    
    if args.output is None:
        plot_outputs_dir = os.path.join(args.out_dir, "plot_outputs")
        os.makedirs(plot_outputs_dir, exist_ok=True)
        args.output = os.path.join(plot_outputs_dir, "sn1_isomerization_comparison.png")
    
    if args.xvg_dirs is not None:
        xvg_dirs = args.xvg_dirs
    elif args.xvg_dir is not None:
        xvg_dirs = [args.xvg_dir]
    else:
        data_dir = os.path.join(args.out_dir, "data")
        xvg_dirs = [
            os.path.join(data_dir, "xvg_files_OFF_HMR"), # temp file name. tbu with l21/c36
            os.path.join(data_dir, "xvg_files_macrog"), 
            os.path.join(data_dir, "xvg_files_slipids")
        ]
    
    valid_dirs = []
    for dir in xvg_dirs:
        if os.path.exists(dir):
            valid_dirs.append(dir)
        else:
            print(f"no xvg.")
    
    xvg_dirs = valid_dirs
    
    if len(xvg_dirs) == 0:
        print("no xvg")
        return 1
    
    labels = args.labels
    if labels is None:
        labels = [os.path.basename(dir).replace('xvg_files_', '') for dir in xvg_dirs]
        if len(labels) < len(xvg_dirs):
            labels.extend([f"Dataset {i+1}" for i in range(len(labels), len(xvg_dirs))])
    
    colors = args.color_scheme
    if len(colors) < len(xvg_dirs):
        default_colors = ['blue', 'green', 'red', 'purple', 'orange', 'brown', 'pink'] # for tests
        colors.extend([default_colors[i % len(default_colors)] for i in range(len(colors), len(xvg_dirs))])
    
    positive_range = tuple(args.positive_range)
    negative_range = tuple(args.negative_range)
    percentages_dict = {dir: [] for dir in xvg_dirs}
    
    for dir in xvg_dirs:
        xvg_files = [os.path.join(dir, f) for f in sorted(os.listdir(dir)) if f.endswith('.xvg')]
        # print(f" {len(xvg_files)} xvg files in {dir}")
        
        for xvg_file in xvg_files:
            combined_percent = calculate_percentage(xvg_file, positive_range, negative_range)
            percentages_dict[dir].append(combined_percent)
            # except Exception as e: # dont really need try
            #     print(f"{xvg_file}: {e}")
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    for dir, color, label in zip(xvg_dirs, colors, labels):
        percentages = percentages_dict[dir]
        if not percentages:
            print(f"No data")
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
    max_bonds = max([len(percentages_dict[dir]) for dir in xvg_dirs]) # xticks=bonds
    ax.set_xticks(np.arange(1, max_bonds + 1))
    fig.tight_layout()
    plt.savefig(args.output, dpi=300)
    print(f"Plot saved to {args.output}")
    
    results_dir = os.path.join(args.out_dir, "results")
    os.makedirs(results_dir, exist_ok=True)
    
    for dir in xvg_dirs:
        percentages = percentages_dict[dir]
        if not percentages:
            continue
            
        output_xvg = os.path.join(results_dir, f'sn1_isomerization_{os.path.basename(dir)}.xvg')
        with open(output_xvg, 'w') as file:
            file.write('# Probability Percentages within Specified Angle Ranges\n')
            file.write(f'# Bond#  Combined({positive_range[0]} to {positive_range[1]} and {negative_range[0]} to {negative_range[1]})\n')
            for i, combined_percent in enumerate(percentages, start=1):
                file.write(f"{i}  {combined_percent:.2f}\n")
        print(f"Results saved to {output_xvg}")
    
    # plt.show()
    
    return 0

def calculate_percentage(xvg_file, positive_range, negative_range):
    aux = mda.auxiliary.XVG.XVGReader(xvg_file)
    data = []
    for step in aux:
        step_data = step.data
        if step_data.ndim == 1:
            step_data = step_data.reshape(-1, 2)
        data.extend(step_data)
    
    df = pd.DataFrame(data, columns=['A', 'B'])
    
    positive_prob = df[(df['A'] >= positive_range[0]) & (df['A'] <= positive_range[1])]['B'].sum()
    negative_prob = df[(df['A'] >= negative_range[0]) & (df['A'] <= negative_range[1])]['B'].sum()
    
    total_prob = df['B'].sum()
    combined_prob = positive_prob + negative_prob
    combined_percentage = (combined_prob / total_prob) * 100 if total_prob != 0 else 0
    
    return combined_percentage

if __name__ == "__main__":
    exit_code = main()
    exit(exit_code)