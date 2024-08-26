import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda

simulation_files = [
    ('Slipids.gro', 'Slipids.xtc'),
    ('2.1_0.9.gro', '2.1_0.9.xtc'),
    ('2.1_1.2.gro', '2.1_1.2.xtc'),
    ('2.0_0.9.gro', '2.0_0.9.xtc'),
    ('2.1_PME.gro', '2.1_PME.xtc'),
]

custom_names = [
    'Slipids',
    '2.1_0.9',
    '2.1_1.2',
    '2.0_0.9',
    '2.1_PME',
]

lipid_resname = 'POPC'

def calculate_area_per_lipid(universe, lipid_resname):
    lipid_group = universe.select_atoms(f"resname {lipid_resname}")
    n_lipids = len(lipid_group.residues)
    
    # if n_lipids == 0:
    #     raise ValueError(f"Resname '{lipid_resname}' is incorrect")
    
    apl = {}

    for ts in universe.trajectory:
        dimensions = ts.dimensions
        apl_frame = (dimensions[0] * dimensions[1] * 2) / n_lipids # APL calc with box area (same as NMR)
        apl[ts.time] = apl_frame

    return apl

def plot_area_per_lipid(simulation_files, lipid_resname, custom_names):
    plt.figure()

    for i, (topology_file, trajectory_file) in enumerate(simulation_files):
        u = mda.Universe(topology_file, trajectory_file)
        
        apl_data = calculate_area_per_lipid(u, lipid_resname)

        # arrays for plotting
        times = np.array(list(apl_data.keys()))
        area_per_lipid = np.array(list(apl_data.values()))

        average_apl = np.mean(area_per_lipid)
        std_apl = np.std(area_per_lipid)

        # plot for data
        plt.plot(times, area_per_lipid, alpha=0.5, label=f"{custom_names[i]}: {average_apl:.2f} ± {std_apl:.2f} Å²")
        
        # overlapping trendline for each
        window_size = 50  # testing idk
        trendline = np.convolve(area_per_lipid, np.ones(window_size)/window_size, mode='valid')
        plt.plot(times[:len(trendline)], trendline, alpha=1.0, color=plt.gca().lines[-1].get_color())

    plt.xlabel('Time (ps)')
    plt.ylabel('Area per lipid (Å²)')
    plt.title('Area per Lipid Over Time')
    plt.xlim(0, 50000)
    plt.legend()
    plt.savefig("50ns_rvdw_APL_trend", dpi=300)
    plt.show()

    # save to file to check with NMR
    with open('area_per_lipid_all.txt', 'w') as f:
        for topology_file, trajectory_file in simulation_files:
            u = mda.Universe(topology_file, trajectory_file)
            apl_data = calculate_area_per_lipid(u, lipid_resname)
            times = np.array(list(apl_data.keys()))
            area_per_lipid = np.array(list(apl_data.values()))
            np.savetxt(f, np.column_stack((times, area_per_lipid)), header=f"{topology_file} Time(ps) Area_per_lipid(Å²)")

plot_area_per_lipid(simulation_files, lipid_resname, custom_names)
