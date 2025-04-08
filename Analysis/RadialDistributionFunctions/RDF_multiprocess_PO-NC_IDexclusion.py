import os
import numpy as np
import MDAnalysis as mda
import pandas as pd
from multiprocessing import Pool
import multiprocessing

index_folder = 'index_files'
xvg_folder = 'xvg_files'

os.makedirs(index_folder, exist_ok=True)
os.makedirs(xvg_folder, exist_ok=True)

# MacRog indices currently
phosphate_indices = 'a P1 a O1 a O2 a O3 a O4'
amine_indices = 'a N1 a C1 a C2 a C3 a C4'

tpr = 'md.tpr'
xtc = 'md.xtc'
amine_name = 'NC'
phosphate_name = 'PO'

density_dict = {}

#range(start,stop,step)
# import parallelTestModule

def f(i):

    index_file_path = os.path.join(index_folder, f'{i}.ndx')

    # build index file
    os.system(f"echo \"ri {i} & {phosphate_indices} \\n name 4 {phosphate_name} \\n ! ri {i} & {amine_indices} \\n name 5 {amine_name} \\n q\" | gmx make_ndx -f {tpr} -o {index_file_path}")
    
    xvg_file_path = os.path.join(xvg_folder, f'{i}.xvg')

    # generate rdf (-xy is 2D rdf, mol_com between two atom groups)
    os.system(f"echo \"\n\" | gmx rdf -f {xtc} -s {tpr} -ref {phosphate_name} -selrpos mol_com -sel {amine_name} -seltype mol_com -n {index_file_path} -xy -o {xvg_file_path}")
    
if __name__ ==  '__main__':
    with Pool(multiprocessing.cpu_count()) as p:
        p.map(f, range(5001, 5129, 1)) # also macrog indices (use 5121,5249 for standard 40w/l for 128 lipids)
 
    for i in range(5001, 5129, 1):

        xvg_file_path = os.path.join(xvg_folder, f'{i}.xvg')
        aux = mda.auxiliary.core.auxreader(xvg_file_path)

        # temp density storage
        column_b = []

        for step in aux:
            step_data = step.data
                
            if step_data.ndim == 1:
                step_data = step_data.reshape(-1, 2)
            for distance, density in step_data:
                if distance not in density_dict:# check for parallel bugs
                    density_dict[distance] = []
                density_dict[distance].append(density)
            # print(density_values_dict[distance])
                
    avg_density_values = {distance: np.mean(densities) for distance, densities in density_dict.items()}
    df = pd.DataFrame(list(avg_density_values.items()), columns=['A', 'B'])

    with open('average_rdf_premade_OFF_2D.xvg', 'w') as f:
        f.write('@TYPE xy\n')
        f.write('@TITLE rdf plot\n')
        f.write('@ xaxis  label "r (nm)"\n')
        f.write('@ yaxis  label "g(r)"\n')
        f.write('@ s0 legend "Data"\n')
        f.write('@ s0 line type 0\n')

        for _, row in df.iterrows():
            f.write(f'{row["A"]} {row["B"]}\n')