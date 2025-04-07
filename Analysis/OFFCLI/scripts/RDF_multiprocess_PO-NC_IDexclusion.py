import os
import numpy as np
import MDAnalysis as mda
import pandas as pd
from multiprocessing import Pool
import multiprocessing
import argparse
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Calculate Radial Distribution Function between phosphate and amine groups on separate lipids")
    parser.add_argument("--tpr", help="Path to TPR file", required=True)
    parser.add_argument("--xtc", help="Path to XTC file", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--output", default="3Drdf_analysis", help="Output file name")
    parser.add_argument("--xy", action="store_true", help="Use 2D RDF calculation (xy-plane only)")
    parser.add_argument("--start-idx", type=int, default=5121, help="Starting residue index + 1")
    parser.add_argument("--end-idx", type=int, default=5249, help="Ending residue index + 1")
    parser.add_argument("--step", type=int, default=1, help="Step size for residue indices (likely not needed)")
    parser.add_argument("--phosphate-indices", default="a P1x a O1x a O2x a O3x a O4x", help="Atom selection for phosphate group")
    parser.add_argument("--amine-indices", default="a N1x a C4x a C5x a C6x a C7x", help="Atom selection for NHG")
    parser.add_argument("--phosphate-name", default="PO", help="Name for phosphate group in index file")
    parser.add_argument("--amine-name", default="NC", help="Name for amine group in index file")
    parser.add_argument("--ndx-group-offset", type=int, default=4, help="Starting group number in index file")
    parser.add_argument("--index-dir")
    parser.add_argument("--xvg-dir")
    return parser.parse_args()

def process_residue(args):
    i, params = args
    out_dir = params["out_dir"]
    tpr = params["tpr"]
    xtc = params["xtc"]
    phosphate_indices = params["phosphate_indices"]
    amine_indices = params["amine_indices"]
    phosphate_name = params["phosphate_name"]
    amine_name = params["amine_name"]
    use_xy = params["use_xy"]
    ndx_group_offset = params["ndx_group_offset"]
    index_dir = params["index_dir"]
    xvg_dir = params["xvg_dir"]

    index_file_path = os.path.join(index_dir, f'{i}.ndx')
    
    po_group = ndx_group_offset
    nc_group = ndx_group_offset + 1
    
    # frequent errors in processing when the selection groups in the calls below are numbered incorrectly
    os.system(f"echo \"ri {i} & {phosphate_indices} \\n name {po_group} {phosphate_name} \\n ! ri {i} & {amine_indices} \\n name {nc_group} {amine_name} \\n q\" | gmx make_ndx -f {tpr} -o {index_file_path}")
    
    xvg_file_path = os.path.join(xvg_dir, f'{i}.xvg')
    
    xy_option = "-xy" if use_xy else ""
    os.system(f"echo \"\n\" | gmx rdf -f {xtc} -s {tpr} -ref {phosphate_name} -selrpos mol_com -sel {amine_name} -seltype mol_com -n {index_file_path} {xy_option} -o {xvg_file_path}")
    
    return xvg_file_path

def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    
    index_folder = args.index_dir if args.index_dir else os.path.join(args.out_dir, "index_files")
    xvg_folder = args.xvg_dir if args.xvg_dir else os.path.join(args.out_dir, "xvg_files")
    os.makedirs(index_folder, exist_ok=True)
    os.makedirs(xvg_folder, exist_ok=True)
    
    params = {
        "out_dir": args.out_dir,
        "tpr": args.tpr,
        "xtc": args.xtc,
        "phosphate_indices": args.phosphate_indices,
        "amine_indices": args.amine_indices,
        "phosphate_name": args.phosphate_name,
        "amine_name": args.amine_name,
        "use_xy": args.xy,
        "ndx_group_offset": args.ndx_group_offset,
        "index_dir": index_folder,
        "xvg_dir": xvg_folder
    }
    
    residue_range = range(args.start_idx, args.end_idx, args.step)

    print(f"Using {'2D (xy-plane)' if args.xy else '3D'} RDF calculation") # maybe make clearer with warn and sleep
    
    with Pool(multiprocessing.cpu_count()) as p:
        p.map(process_residue, [(i, params) for i in residue_range])
    
    density_dict = {}
    
    for i in residue_range:
        xvg_file_path = os.path.join(xvg_folder, f'{i}.xvg')
        
        try:
            aux = mda.auxiliary.core.auxreader(xvg_file_path)
            
            for step in aux:
                step_data = step.data
                
                if step_data.ndim == 1:
                    step_data = step_data.reshape(-1, 2)
                for distance, density in step_data:
                    if distance not in density_dict:
                        density_dict[distance] = []
                    density_dict[distance].append(density)
        except Exception as e:
            print(f"error: {xvg_file_path} / {e}")
    
    avg_density_values = {distance: np.mean(densities) for distance, densities in density_dict.items()}
    df = pd.DataFrame(list(avg_density_values.items()), columns=['Distance', 'Density'])
    df = df.sort_values('Distance')
    csv_output = os.path.join(args.out_dir, f"{args.output}.csv")
    df.to_csv(csv_output, index=False)
    xvg_output = os.path.join(args.out_dir, f"{args.output}.xvg")
    
    with open(xvg_output, 'w') as f: # write xvg
        f.write('@TYPE xy\n')
        f.write(f'@TITLE RDF between {args.phosphate_name} and {args.amine_name}\n')
        f.write(f'@ xaxis label "r (nm)"\n')
        f.write('@ yaxis label "g(r)"\n')
        f.write('@ s0 legend "Data"\n')
        f.write('@ s0 line type 0\n')
        
        for _, row in df.iterrows():
            f.write(f'{row["Distance"]} {row["Density"]}\n')
    
    print(f"Saved to {xvg_output}")

if __name__ == '__main__':
    main()