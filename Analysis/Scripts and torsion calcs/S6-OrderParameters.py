import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import time

trajectory = md.load('md.xtc', top='md.gro')

acyl_chain_patterns = [
    # Openff
    ['C8x', 'C9x', 'C10x', 'C11x', 'C12x', 'C13x', 'C14x', 'C15x', 'C16x', 'C17x', 'C18x', 'C19x', 'C20x', 'C21x', 'C22x', 'C23x', 'C24x', 'C25x'],
    ['C27x', 'C28x', 'C29x', 'C30x', 'C31x', 'C32x', 'C33x', 'C34x', 'C35x', 'C36x', 'C37x', 'C38x', 'C39x', 'C40x', 'C41x', 'C42x']
    # MacRog
    # ['C8', 'C9', 'C13', 'C14', 'C15', 'C16', 'C17', 'C18', 'C19', 'C20', 'C21', 'C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28'],
    # ['C11', 'C12', 'C29', 'C30', 'C31', 'C32', 'C33', 'C34', 'C35', 'C36', 'C37', 'C38', 'C39', 'C40', 'C41', 'C42']
    # Slipids (HAS 122 LIPIDS)
    # ['C21', 'C22', 'C23', 'C24', 'C25', 'C26', 'C27', 'C28', 'C29', 'C210', 'C211', 'C212', 'C213', 'C214', 'C215', 'C216', 'C217', 'C218'],
    # ['C31', 'C32', 'C33', 'C34', 'C35', 'C36', 'C37', 'C38', 'C39', 'C310', 'C311', 'C312', 'C313', 'C314', 'C315', 'C316']
]

def separate_acyl_chains(topology, patterns):
    chains = []                             # list of lists
    for residue in topology.residues:
        for pattern in patterns:
            chain = [atom.index for atom in residue.atoms if atom.name in pattern] # find residue index that matches a pattern
            if len(chain) == len(pattern):  # full pattern match for single acyl chain
                chains.append(chain)
    return chains

def compute_average_positions(positions, chains):
    average_positions = [] #pos array
    for chain in chains:
        chain_positions = positions[chain, :2] # chain index to get atoms in position array and just x y
        avg_position = np.mean(chain_positions, axis=0) # axis=0 averages the x and y in each ROW
        average_positions.append(avg_position)
    return np.array(average_positions)

# get positions of the final frame + calls
final_positions = trajectory.xyz[-1]
acyl_chains = separate_acyl_chains(trajectory.topology, acyl_chain_patterns)
average_positions = compute_average_positions(final_positions, acyl_chains)
# ensure average_positions is a 2D array with shape (num_atoms, 2)      (bug check)
average_positions_2d = average_positions.reshape(-1, 2)

# DEBUG
print(f"Shape of average_positions: {average_positions.shape}")       # 256
print(f"Shape of average_positions_2d: {average_positions_2d.shape}") # 256
unique_positions = np.unique(average_positions_2d, axis=0)            # 256
print(f"Number of unique positions: {unique_positions.shape[0]}")     # 256
print(f"Number of acyl chains: {len(acyl_chains)}")                   # 256
for i, chain in enumerate(acyl_chains[:3]):     # only 3
    print(f"Acyl Chain {i}: Indices: {chain}")  # 16 or 18 indices

# calculate the hexagonal order parameter (S6) for a given index
def calculate_S6_for_index(i, positions, cutoff=0.65):

    # start_time = time.time()

    num_atoms = positions.shape[0]
    # print(f"Total atoms for loop: {num_atoms}") # 256

    neighbors = []
    for j in range(num_atoms):
        if i != j:
            dx = abs(positions[i, 0] - positions[j, 0])
            dy = abs(positions[i, 1] - positions[j, 1])
            if dx >= cutoff or dy >= cutoff: # optimization: skip most outside cutoff
                continue
            distance = np.sqrt(dx**2 + dy**2) # vector length
            if distance < cutoff:
                neighbors.append(j)
    
    if len(neighbors) == 0:
        return 0.0
    
    psi_jk = [] # list of angles between i and all js for sum
    for neighbor in neighbors:
        delta = positions[neighbor] - positions[i]  # get vector
        angle = np.arctan2(delta[1], delta[0]) # find angle (radians???) from vector (between x-axis and vector with origin as reference)
        psi_jk.append(angle)
    
    psi_jk = np.array(psi_jk) # convert list to array
    S6 = (1/6) * np.sqrt((np.sum(np.cos(6 * psi_jk)))**2 + (np.sum(np.sin(6 * psi_jk)))**2)

    # print(f"Total execution time: --- { (time.time() - start_time):.2f} seconds ---\n\n" )

    return S6

def calculate_S6_multiprocessing(positions, cutoff=0.65):

    num_atoms = positions.shape[0]
    # q = [(i, positions, cutoff) for i in range(num_atoms)]
    # print(f"Q list: {len(q)}")

    with Pool(4) as p:
        S6 = list(tqdm(p.starmap(calculate_S6_for_index, [(i, positions, cutoff) for i in range(num_atoms)]), total=num_atoms)) # total=num_atoms for progress bar
    return np.array(S6)

if __name__ == '__main__':
    S6_last_frame = calculate_S6_multiprocessing(average_positions_2d) # main call with averages

    # outside 0 to 1 (DEBUG)
    out_of_range_S6 = S6_last_frame[(S6_last_frame < 0) | (S6_last_frame > 1)]
    if len(out_of_range_S6) > 0:
        print(f"S6 values out of range: {out_of_range_S6}")


    average_S6 = np.mean(S6_last_frame)
    print(f"Average S6 order parameter: {average_S6}")
    std_S6 = np.std(S6_last_frame)
    print(f"Standard deviation of S6 order parameter: {std_S6}")
    sem_S6 = std_S6/np.sqrt(S6_last_frame.size)
    print(f"Standard error of S6 order parameter: {sem_S6}")
    print(S6_last_frame.size)

    def plot_hexagonal_order_parameter(positions, S6, title):

        plt.figure()
        
        # only x and y coordinates for Voronoi diagram (voronoi cell is closer to point than anything else)
        positions_2d = positions[:, :2] 
        
        vor = Voronoi(positions_2d)
        # voronoi_plot_2d(vor, show_vertices=False, line_colors='k', line_width=0.5, line_alpha=0.6, point_size=2)
        for simplex in vor.ridge_vertices: # loop each boundary
            simplex = np.asarray(simplex) # index list to array
            if np.all(simplex >= 0): # check non-negative (within plot boundary)
                plt.plot(vor.vertices[simplex, 0], vor.vertices[simplex, 1], 'k-', lw=0.5, alpha=0.6) # x then y, k- = black

        #average pos for point, s6 color
        scatter = plt.scatter(positions_2d[:, 0], positions_2d[:, 1], c=S6, cmap='coolwarm_r', s=15, edgecolors=None) # "_r" flips color bar & c=S6 should scale colors...
        plt.colorbar(scatter, label='$S_6$', extend='both')
        
        plt.xlim(np.min(positions_2d[:, 0]), np.max(positions_2d[:, 0])) # x range
        plt.ylim(np.min(positions_2d[:, 1]), np.max(positions_2d[:, 1])) # y range
        plt.gca().set_aspect('equal', 'box') # equal xy scaling
        plt.title(title)

        plt.savefig(title, dpi=300)
        plt.show()

    top_positions = average_positions_2d[:128] # CHANGE FOR SLIPIDS TO 122
    top_S6 = S6_last_frame[:128]
    plot_hexagonal_order_parameter(top_positions, top_S6, 'OpenFF (320K) Top Monolayer Hexagonal Order Parameter ($S_6$)')

    bottom_positions = average_positions_2d[128:]
    bottom_S6 = S6_last_frame[128:]
    plot_hexagonal_order_parameter(bottom_positions, bottom_S6, 'OpenFF (320K) Bottom Monolayer Hexagonal Order Parameter ($S_6$)')