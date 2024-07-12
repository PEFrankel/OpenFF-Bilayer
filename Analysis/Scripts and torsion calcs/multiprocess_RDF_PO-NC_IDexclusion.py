import os
import numpy as np
import MDAnalysis as mda
import pandas as pd
from multiprocessing import Pool
import multiprocessing

# NEED TO CHANGE INDICES, FILE NAMES, & INDEX # IN ECHO CALLS
# NEED TO CHANGE INDICES, FILE NAMES, & INDEX # IN ECHO CALLS
# NEED TO CHANGE INDICES, FILE NAMES, & INDEX # IN ECHO CALLS
# NEED TO CHANGE INDICES, FILE NAMES, & INDEX # IN ECHO CALLS

# organization
index_folder = 'index_files'
xvg_folder = 'xvg_files'

os.makedirs(index_folder, exist_ok=True)
os.makedirs(xvg_folder, exist_ok=True)

# MacRog indices currently
phosphate_indices = 'a P1 a O1 a O2 a O3 a O4'
amine_indices = 'a N1 a C1 a C2 a C3 a C4'

tpr = 'md.tpr'
xtc = 'md.xtc'      # assumes tpr and trj have these exact names (functional for my purposes)
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
        p.map(f, range(5001, 5129, 1))
 
    for i in range(5001, 5129, 1):

        xvg_file_path = os.path.join(xvg_folder, f'{i}.xvg')

        # read the XVG using mda
        aux = mda.auxiliary.core.auxreader(xvg_file_path)

        # temp density storage
        column_b = []

        # iterate for each XVG row from mda
        for step in aux:
            step_data = step.data                   # property 'data' of 'XVGStep' object has no setter
                
            if step_data.ndim == 1:  # check step_data is paired
                step_data = step_data.reshape(-1, 2)

            # iterate over each row in step_data
            for distance, density in step_data:
                if distance not in density_dict: # check if the distance is unique
                    density_dict[distance] = []  # if yes, initialize a list to store density values for that distance
                density_dict[distance].append(density) # add the density value to the list for that distance
            # print(density_values_dict[distance])
                
    # average column b at each row for all files
    avg_density_values = {distance: np.mean(densities) for distance, densities in density_dict.items()}

    # build a dataframe with distance and averaged densities
    df = pd.DataFrame(list(avg_density_values.items()), columns=['A', 'B'])

    # df.to_csv('AVERAGE_RDF.csv', index=False)

    with open('average_rdf_premade_OFF_2D.xvg', 'w') as f:
        f.write('@TYPE xy\n')  # Specifies the data type
        f.write('@TITLE rdf plot\n')  # Specifies the title
        f.write('@ xaxis  label "r (nm)"\n')  # Specifies the label for the x-axis
        f.write('@ yaxis  label "g(r)"\n')  # Specifies the label for the y-axis
        f.write('@ s0 legend "Data"\n')  # Specifies the legend
        f.write('@ s0 line type 0\n')  # Specifies the line type

        # Write each row of the DataFrame to the XVG file
        for _, row in df.iterrows():
            f.write(f'{row["A"]} {row["B"]}\n')