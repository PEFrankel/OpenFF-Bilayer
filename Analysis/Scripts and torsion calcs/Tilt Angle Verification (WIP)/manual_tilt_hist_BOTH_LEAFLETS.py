import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import get_context, cpu_count # :)
from scipy import stats






#  ADD COMMENTS     





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

def calculate_normal_vector(residues, phosphorus_atom='P'):
    vectors = []
    for res in residues:
        phosphorus = res.atoms.select_atoms(f"name {phosphorus_atom}")
        if len(phosphorus) > 0:
            vectors.append(phosphorus.positions[0])
    # if len(vectors) == 0:
    #     print("No vectors found for any simulations")                          # no bug -- check this
    #     return np.array([0, 0, 1])
    avg_position = np.mean(vectors, axis=0)
    normal_vector = np.mean(vectors - avg_position, axis=0)
    normal_vector /= np.linalg.norm(normal_vector)
    return normal_vector

def calculate_tilt_angle(residue, normal_vector, phosphorus_atom='P', nitrogen_atom='N'):
    phosphorus = residue.atoms.select_atoms(f"name {phosphorus_atom}")
    nitrogen = residue.atoms.select_atoms(f"name {nitrogen_atom}")
    if len(phosphorus) == 0 or len(nitrogen) == 0:
        return None
    vector = nitrogen.positions[0] - phosphorus.positions[0]
    cosine_angle = np.dot(vector, normal_vector) / (np.linalg.norm(vector) * np.linalg.norm(normal_vector))
    tilt_angle = np.arccos(cosine_angle) * (180.0 / np.pi)
    return tilt_angle

def process_frame(ts, leaflet1, leaflet2, phosphorus_atom='P', nitrogen_atom='N'):
    normal_vector_leaflet1 = calculate_normal_vector(leaflet1, phosphorus_atom)
    normal_vector_leaflet2 = calculate_normal_vector(leaflet2, phosphorus_atom)
    tilt_angles_leaflet1_frame = [calculate_tilt_angle(res, normal_vector_leaflet1, phosphorus_atom, nitrogen_atom) for res in leaflet1]
    tilt_angles_leaflet2_frame = [calculate_tilt_angle(res, normal_vector_leaflet2, phosphorus_atom, nitrogen_atom) for res in leaflet2]
    tilt_angles_leaflet1_frame = [angle for angle in tilt_angles_leaflet1_frame if angle is not None]
    tilt_angles_leaflet2_frame = [angle for angle in tilt_angles_leaflet2_frame if angle is not None]
    return tilt_angles_leaflet1_frame, tilt_angles_leaflet2_frame

def process_simulation_wrapper(args):
    topology, trajectory, phosphorus_atom, nitrogen_atom = args
    context = get_context("spawn")                                      # TERRIBLE -- FIX
    universe = mda.Universe(topology, trajectory)
    try:
        lipid_selection = universe.select_atoms(f'resname POPC')
        leaflet1, leaflet2 = assign_leaflets(lipid_selection)
        with context.Pool(processes=cpu_count()) as pool:
            results = pool.starmap(process_frame, [(ts, leaflet1, leaflet2, phosphorus_atom, nitrogen_atom) for ts in universe.trajectory])
    finally:
        pass  # with covers this
    
    tilt_angles_leaflet1 = [angle for frame in results for angle in frame[0]]
    tilt_angles_leaflet2 = [angle for frame in results for angle in frame[1]]
    avg_tilt_angle_leaflet1 = np.mean(tilt_angles_leaflet1)
    sem_tilt_angle_leaflet1 = stats.sem(tilt_angles_leaflet1, nan_policy='omit') if tilt_angles_leaflet1 else np.nan
    avg_tilt_angle_leaflet2 = np.mean(tilt_angles_leaflet2)
    sem_tilt_angle_leaflet2 = stats.sem(tilt_angles_leaflet2, nan_policy='omit') if tilt_angles_leaflet2 else np.nan
    return avg_tilt_angle_leaflet1, sem_tilt_angle_leaflet1, avg_tilt_angle_leaflet2, sem_tilt_angle_leaflet2

if __name__ == '__main__':
    simulations = [
        ('md_OFF.tpr', 'md_OFF.xtc', 'P1x', 'N1x'),
        ('md_macrog.tpr', 'md_macrog.xtc', 'P1', 'N1'),
        ('md_slipids.tpr', 'md_slipids.xtc', 'P', 'N')
    ]

    context = get_context("spawn")  # :(
    with context.Pool(processes=cpu_count()) as pool:
        results_parallel = pool.map(process_simulation_wrapper, simulations)

    for i, (topology, trajectory, phosphorus_atom, nitrogen_atom) in enumerate(simulations):
        sim_name = f"Simulation {i+1}"
        avg_tilt_angle_leaflet1, sem_tilt_angle_leaflet1, avg_tilt_angle_leaflet2, sem_tilt_angle_leaflet2 = results_parallel[i]
        print(f"{sim_name}")
        print(f"  Average tilt angle for leaflet 1: {avg_tilt_angle_leaflet1:.2f} degrees ± {sem_tilt_angle_leaflet1:.2f}")
        print(f"  Average tilt angle for leaflet 2: {avg_tilt_angle_leaflet2:.2f} degrees ± {sem_tilt_angle_leaflet2:.2f}")
        print()

    # plots + main
    fig, ax = plt.subplots()

    bins = np.linspace(0, 180, 36)

    for i, (topology, trajectory, phosphorus_atom, nitrogen_atom) in enumerate(simulations):
        universe = mda.Universe(topology, trajectory)
        lipid_selection = universe.select_atoms(f'resname POPC')
        leaflet1, leaflet2 = assign_leaflets(lipid_selection)
        tilt_angles_leaflet1 = [angle for frame in results_parallel[i] for angle in frame[0]]
        tilt_angles_leaflet2 = [angle for frame in results_parallel[i] for angle in frame[1]]
        sim_name = f"Simulation {i+1}"
        ax.hist(tilt_angles_leaflet1, bins=bins, density=True, alpha=0.6, color='blue', label=f"{sim_name} Leaflet 1", histtype='stepfilled', linewidth=2)
        ax.hist(tilt_angles_leaflet2, bins=bins, density=True, alpha=0.6, color='orange', label=f"{sim_name} Leaflet 2", histtype='step', linestyle='--', linewidth=2)

    ax.set_title("Tilt Angle of the PN Vector")
    ax.set_xlabel("Degrees")
    ax.set_ylabel("Probability")
    ax.legend()

    plt.show()
