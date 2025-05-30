import os
import numpy as np
import MDAnalysis as mda
import pandas as pd
from collections import defaultdict 

index_folder = 'index_files'
xvg_folder = 'xvg_files'

os.makedirs(index_folder, exist_ok=True)
os.makedirs(xvg_folder, exist_ok=True)

# OFF indices currently
phosphate_indices = 'a P1x a O1x a O2x a O3x a O4x'
amine_indices = 'a N1x a C4x a C5x a C6x a C7x'

tpr = 'md.tpr'
xtc = 'md.xtc'
amine_name = 'NC'
phosphate_name = 'PO'

density_dict = defaultdict(list)

# range(start, stop, step)
for i in range(5121, 5249):  # +1 to stop; stops one prior: / for each resid (128 in this case)
    index_file_path = os.path.join(index_folder, f'{i}.ndx')
    xvg_file_path = os.path.join(xvg_folder, f'{i}.xvg')
    
    if not os.path.exists(xvg_file_path):
        print(f"Skipping {xvg_file_path}: File does not exist")
        continue
    
    try:
        aux = mda.auxiliary.XVG.XVGReader(xvg_file_path)
    except Exception as e:
        print(f"Error reading {xvg_file_path}: {e}")
        continue

    for step in aux:
        step_data = step.data

        if step_data.ndim == 1:
            step_data = step_data.reshape(-1, 2)
        for distance, density in step_data:
            density_dict[distance].append(density)

avg_density_values = {distance: np.mean(densities) for distance, densities in density_dict.items()}
df = pd.DataFrame(list(avg_density_values.items()), columns=['A', 'B'])

with open('average_rdf_OFF.xvg', 'w') as f:
    f.write('@TYPE xy\n')
    f.write('@TITLE rdf plot\n')
    f.write('@ xaxis  label "r (nm)"\n')
    f.write('@ yaxis  label "g(r)"\n')  
    f.write('@ s0 legend "Data"\n')
    f.write('@ s0 line type 0\n')

    for _, row in df.iterrows():
        f.write(f'{row["A"]} {row["B"]}\n')

print("Averaged RDF data written to average_rdf_OFF.xvg")
