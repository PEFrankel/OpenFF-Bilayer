import os
import numpy as np
import MDAnalysis as mda
import pandas as pd
from multiprocessing import Pool
import multiprocessing
import argparse
import sys

# For bugs ensure changes were properly made to: 
# 1. Loop and p.map indices: 'for i in range(5001, 5129, 1):' to account for your system range.
# 2. Ensure your input and output files are named correctly.
# 3. Index numbers in echo calls: i.e. '\\n name 4' where 4 is the next index for gmx ndx, followed by 5.
#                                 Your system may have more or less default groups in gmx ndx depending or
#                                 FF and/or composition. See 'ndx_group_offset' if adjusting this parameter.

def parse_args():
    parser = argparse.ArgumentParser(description="Calculate Radial Distribution Function (RDF) between phosphate and amine groups on separate lipids")
    parser.add_argument("--tpr", help="Path to TPR topology file")
    parser.add_argument("--xtc", help="Path to XTC trajectory file")
    parser.add_argument("--output", default="3Drdf_analysis", help="Output file name (without extension)")
    parser.add_argument("--xy", action="store_true", help="Use 2D RDF calculation (xy-plane only)")
    parser.add_argument("--start-idx", type=int, default=5001, help="Starting residue index")
    parser.add_argument("--end-idx", type=int, default=5129, help="Ending residue index (exclusive)")
    parser.add_argument("--step", type=int, default=1, help="Step size for residue indices (likely not needed)")
    parser.add_argument("--phosphate-indices", default="a P1x a O1x a O2x a O3x a O4x", help="Atom selection for phosphate group within gmx ndx")
    parser.add_argument("--amine-indices", default="a N1x a C4x a C5x a C6x a C7x", help="Atom selection for amine group")
    parser.add_argument("--phosphate-name", default="PO", help="Name for phosphate group in index file")
    parser.add_argument("--amine-name", default="NC", help="Name for amine group in index file")
    parser.add_argument("--ndx-group-offset", type=int, default=4, help="Starting group number in index file")
    return parser.parse_args()

def process_residue(args):
    i, params = args
    index_folder = params["index_folder"]
    xvg_folder = params["xvg_folder"]
    tpr = params["tpr"]
    xtc = params["xtc"]
    phosphate_indices = params["phosphate_indices"]
    amine_indices = params["amine_indices"]
    phosphate_name = params["phosphate_name"]
    amine_name = params["amine_name"]
    use_xy = params["use_xy"]
    ndx_group_offset = params["ndx_group_offset"]

    index_file_path = os.path.join(index_folder, f'{i}.ndx')
    
    # Build index file with group numbers based on the offset
    po_group = ndx_group_offset
    nc_group = ndx_group_offset + 1
    
    # Create index file command
    os.system(f"echo \"ri {i} & {phosphate_indices} \\n name {po_group} {phosphate_name} \\n ! ri {i} & {amine_indices} \\n name {nc_group} {amine_name} \\n q\" | gmx make_ndx -f {tpr} -o {index_file_path}")
    
    xvg_file_path = os.path.join(xvg_folder, f'{i}.xvg')
    
    # Generate RDF command
    xy_option = "-xy" if use_xy else ""
    os.system(f"echo \"\n\" | gmx rdf -f {xtc} -s {tpr} -ref {phosphate_name} -selrpos mol_com -sel {amine_name} -seltype mol_com -n {index_file_path} {xy_option} -o {xvg_file_path}")
    
    return xvg_file_path

def main():
    args = parse_args()
    
    # Get parameters from arguments or environment variables
    tpr = args.tpr or os.environ.get("TPR_FILE", "md.tpr")
    xtc = args.xtc or os.environ.get("XTC_FILE", "md.xtc")
    output_name = args.output
    use_xy = args.xy
    start_idx = args.start_idx
    end_idx = args.end_idx
    step = args.step
    phosphate_indices = args.phosphate_indices
    amine_indices = args.amine_indices
    phosphate_name = args.phosphate_name
    amine_name = args.amine_name
    ndx_group_offset = args.ndx_group_offset
    
    # Check if input files exist
    if not os.path.exists(tpr):
        print(f"Error: TPR file '{tpr}' does not exist")
        sys.exit(1)
    if not os.path.exists(xtc):
        print(f"Error: XTC file '{xtc}' does not exist")
        sys.exit(1)
    
    # Create directories for output
    index_folder = 'index_files'
    xvg_folder = 'xvg_files'
    os.makedirs(index_folder, exist_ok=True)
    os.makedirs(xvg_folder, exist_ok=True)
    
    # Parameters to pass to worker function
    params = {
        "index_folder": index_folder,
        "xvg_folder": xvg_folder,
        "tpr": tpr,
        "xtc": xtc,
        "phosphate_indices": phosphate_indices,
        "amine_indices": amine_indices,
        "phosphate_name": phosphate_name,
        "amine_name": amine_name,
        "use_xy": use_xy,
        "ndx_group_offset": ndx_group_offset
    }
    
    # Create residue range
    residue_range = range(start_idx, end_idx, step)
    print(f"Processing residue indices from {start_idx} to {end_idx-1} with step {step}")
    print(f"Using {'2D (xy-plane)' if use_xy else '3D'} RDF calculation")
    
    # Parallelization
    with Pool(multiprocessing.cpu_count()) as p:
        p.map(process_residue, [(i, params) for i in residue_range]) # (ensure residue_range is consistent with your system)
    
    # Dictionary to store density values for each distance
    density_dict = {}
    
    # Process each XVG file
    for i in residue_range: # (ensure residue_range is consistent with your system)
        xvg_file_path = os.path.join(xvg_folder, f'{i}.xvg')
        
        try:
            aux = mda.auxiliary.core.auxreader(xvg_file_path)
            
            # Iterate for each XVG row
            for step in aux:
                step_data = step.data
                
                if step_data.ndim == 1:
                    step_data = step_data.reshape(-1, 2)
                
                # Iterate over each row in step_data
                for distance, density in step_data:
                    if distance not in density_dict:
                        density_dict[distance] = []
                    density_dict[distance].append(density)
        except Exception as e:
            print(f"Warning: Error processing {xvg_file_path}: {e}")
    
    # Calculate average density values
    avg_density_values = {distance: np.mean(densities) for distance, densities in density_dict.items()}
    df = pd.DataFrame(list(avg_density_values.items()), columns=['Distance', 'Density'])
    df = df.sort_values('Distance')
    
    # Save to CSV (if interested)
    csv_output = f"{output_name}.csv"
    df.to_csv(csv_output, index=False)
    print(f"Saved CSV data to {csv_output}")
    
    # Save to XVG for visualization
    xvg_output = f"{output_name}.xvg"
    dimension_label = "r_xy (nm)" if use_xy else "r (nm)"
    
    with open(xvg_output, 'w') as f:
        f.write('@TYPE xy\n')
        f.write(f'@TITLE RDF between {phosphate_name} and {amine_name}\n')
        f.write(f'@ xaxis label "{dimension_label}"\n')
        f.write('@ yaxis label "g(r)"\n')
        f.write('@ s0 legend "Data"\n')
        f.write('@ s0 line type 0\n')
        
        for _, row in df.iterrows():
            f.write(f'{row["Distance"]} {row["Density"]}\n')
    
    print(f"Saved XVG plot data to {xvg_output}")
    print(f"Analysis complete!")

if __name__ == '__main__':
    main()