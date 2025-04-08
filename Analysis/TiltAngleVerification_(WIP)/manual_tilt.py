import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

def assign_leaflets(atom_group):

    leaflet1 = []
    leaflet2 = []

    com_z_values = [res.atoms.center_of_mass()[2] for res in atom_group.residues]
    z_cutoff = np.mean(com_z_values)

    for res, com_z in zip(atom_group.residues, com_z_values):
        if com_z < z_cutoff:
            leaflet1.append(res)
        else:
            leaflet2.append(res)
    
    return leaflet1, leaflet2

def calculate_tilt_angle(residue, phosphorus_atom='P', nitrogen_atom='N'):

    phosphorus = residue.atoms.select_atoms(f"name {phosphorus_atom}")
    nitrogen = residue.atoms.select_atoms(f"name {nitrogen_atom}")

    if len(phosphorus) == 0 or len(nitrogen) == 0:
        return None

    vector = nitrogen.positions[0] - phosphorus.positions[0]
    z_vector = np.array([0, 0, 1])
    cosine_angle = np.dot(vector, z_vector) / (np.linalg.norm(vector) * np.linalg.norm(z_vector))
    tilt_angle = np.arccos(cosine_angle) * (180.0 / np.pi)

    return tilt_angle

def read_sim(topology, trajectory, phosphorus_atom='P', nitrogen_atom='N'): # for trj read/avg

    universe = mda.Universe(topology, trajectory)
    lipid_selection = universe.select_atoms(f'resname POPC')
    leaflet1, leaflet2 = assign_leaflets(lipid_selection)

    tilt_angles_leaflet1 = [calculate_tilt_angle(res, phosphorus_atom, nitrogen_atom) for res in leaflet1]
    tilt_angles_leaflet2 = [calculate_tilt_angle(res, phosphorus_atom, nitrogen_atom) for res in leaflet2]
    tilt_angles_leaflet1 = [angle for angle in tilt_angles_leaflet1 if angle is not None]
    tilt_angles_leaflet2 = [angle for angle in tilt_angles_leaflet2 if angle is not None]

    avg_tilt_angle_leaflet1 = np.mean(tilt_angles_leaflet1) if tilt_angles_leaflet1 else None
    avg_tilt_angle_leaflet2 = np.mean(tilt_angles_leaflet2) if tilt_angles_leaflet2 else None

    return avg_tilt_angle_leaflet1, avg_tilt_angle_leaflet2

# simulation files and paths -- adjust as necessary
simulations = [
    ('md_OFF.tpr', 'md_OFF.xtc', 'P1x', 'N1x', 'OpenFF', '#0333b0'),
    ('md_macrog.tpr', 'md_macrog.xtc', 'P1', 'N1', 'MacRog', 'green'),
    ('md_slipids.tpr', 'md_slipids.xtc', 'P', 'N', 'Slipids', '#ee7f17')
]

results = []
for i, (topology, trajectory, phosphorus_atom, nitrogen_atom, sim_name, color) in enumerate(simulations):
    avg_tilt_angle_leaflet1, avg_tilt_angle_leaflet2 = read_sim(topology, trajectory, phosphorus_atom=phosphorus_atom, nitrogen_atom=nitrogen_atom)
    results.append((sim_name, color, avg_tilt_angle_leaflet1, avg_tilt_angle_leaflet2))

fig, ax = plt.subplots()

for sim_name, color, avg_tilt_angle_leaflet1, avg_tilt_angle_leaflet2 in results:
    if avg_tilt_angle_leaflet1 is not None:
        ax.bar(f"{sim_name} Leaflet 1", avg_tilt_angle_leaflet1, color=color, alpha=0.6, label=f"{sim_name} Leaflet 1")
    if avg_tilt_angle_leaflet2 is not None:
        ax.bar(f"{sim_name} Leaflet 2", avg_tilt_angle_leaflet2, color=color, alpha=0.9, label=f"{sim_name} Leaflet 2")

ax.set_title("Tilt Angle of the PN Vector")
ax.set_ylabel("Degrees")
ax.legend()

plt.show()