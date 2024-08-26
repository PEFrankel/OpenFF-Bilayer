import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import time
import math
from functools import reduce

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

acyl_chains = separate_acyl_chains(trajectory.topology, acyl_chain_patterns)

def compute_average_positions(positions, chains):
    # start_time = time.time()
    average_positions = []
    for chain in chains:
        chain_positions = positions[:, chain, :2] # half frames of chain index to get atoms in position array and just x y
        avg_position = np.mean(chain_positions, axis=1) # axis=0 averages the x and y in each ROW
        average_positions.append(avg_position)
    # print(f"Total execution time: --- { (time.time() - start_time):.2f} seconds ---\n\n" )
    return np.array(average_positions)

# for plot
final_frame_positions = trajectory.xyz[-1]
final_positions = compute_average_positions(final_frame_positions[np.newaxis, :], acyl_chains)
final_positions_2d = final_positions.reshape(-1, 2)

# get positions of last half of frames + calls
num_frames = trajectory.n_frames
last_half_frames = trajectory.xyz[(31 * num_frames) // 32:] # gets all positions (2501, 32512, 3) for the last 1/8 -> '[num_frames // 2:]' for last half
print(last_half_frames.size)
print(last_half_frames.shape)
print("^^^^^^^^^^^^^^")

average_positions = compute_average_positions(last_half_frames, acyl_chains)
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


def find_neighbors(i, positions, cutoff=0.65): # BROADCASTING OPTIMIZATION (get rid of loop and less distance calculations)
    # start_time = time.time()
    # (
    delta = np.abs(positions - positions[i]) # vectorized: check distance between i and all atoms (only access the lists twice and with less specificity)
    
    within_cutoff = np.all(delta < cutoff, axis=1) # check cut off for x or y (np.all gives boolean array)
    within_cutoff[i] = False # i != j
    # ) avoids loop through atoms

    distances_squared = np.sum(delta[within_cutoff]**2, axis=1) # get squared distance for all potential neighbors to avoid loop (both x,y like before)
    neighbors = np.where(np.sqrt(distances_squared) < cutoff)[0] # check length again and return chain indices relative to those within cut off
    
    neighbor_indices = np.where(within_cutoff)[0][neighbors] # checks chains that pass the initial, then final cut off tests/indices need to be corrected to refer to the original array?

    # Vectorized operations: delta array is computed for all atoms relative to the atom at index i, 
    # then apply the cutoff to filter out atoms that are outside the cut off 
    # This avoids the need to loop over all atoms explicitly

    # Broadcasting: delta array is used to calculate squared distances for all potential neighbors in a single operation, 
    # avoiding loops through each atom / the smaller array is “broadcast” across the larger array so that they have compatible shapes

    # filter neighbors based on the actual Euclidean distance only after applying the simpler x and y cut offs
    # print(f"Total execution time: --- neighbors size = {len(neighbor_indices):,} ----{ (time.time() - start_time):.2f} seconds ---\n\n" )


    # print(f"Atom {i}: Neighbors found: {neighbor_indices.size}")
    # print(f"Neighbor indices: {neighbor_indices}")
    # print(f"Delta values: {delta}")
    # print(f"Distances squared: {distances_squared}")
    # print(f"Distances: {np.sqrt(distances_squared)}")

    return neighbor_indices


# calculate the hexagonal order parameter (S6) for a given index
def calculate_S6_for_index(i, positions, cutoff=0.65):
    neighbors = find_neighbors(i, positions, cutoff)
    
    if len(neighbors) == 0:
        return 0.0
    
    # delta = positions[neighbors] - positions[i]
    # angles = np.arctan2(delta[:, 1], delta[:, 0]) # delta is array of 2d vectors between i and neighbors (num_neighbors, 2)

    # cos_sum = np.sum(np.cos(6 * angles))
    # sin_sum = np.sum(np.sin(6 * angles))

    # S6 = (1/6) * np.sqrt(cos_sum**2 + sin_sum**2)

    ######################################################################################
    ######################################################################################
    S6_values = []

    # Calculate S6 for each neighbor
    for neighbor in neighbors:
        delta = positions[neighbor] - positions[i]
        angle = np.arctan2(delta[1], delta[0])
        # print(f"Index {i}, Neighbor {neighbor}: Angle: {angle}")
        # Compute the S6 value for this particular neighbor and add it to the list
        S6 = (1/6) * np.sqrt((np.sum(np.cos(6 * angle)))**2 + (np.sum(np.sin(6 * angle)))**2)
        S6_values.append(S6)
        # print(f"Index {i}, Neighbor {neighbor}: S6: {S6}")
    # Average the S6 values for all neighbors
    average_S6 = np.mean(S6_values)

    ######################################################################################
    ######################################################################################

    # psi_jk = [] # list of angles between i and all js for sum
    # for neighbor in neighbors:
    #     delta = positions[neighbor] - positions[i]  # get displacement vector
    #     angle = np.arctan2(delta[1], delta[0]) # find angle (radians???) from vector (between x-axis and vector with origin as reference)
    #     psi_jk.append(angle)
    
    # psi_jk = np.array(psi_jk) # convert list to array
    # S6 = (1/6) * np.sqrt((np.sum(np.cos(6 * psi_jk)))**2 + (np.sum(np.sin(6 * psi_jk)))**2)

    return average_S6


def process_index(i, positions, cutoff):
    S6 = calculate_S6_for_index(i, positions, cutoff)
    if i % 25000 == 0 and i != 0:
        print(f"Processed {i} chains out of {positions.shape[0]}")
    return S6

def calculate_S6_multiprocessing(positions, cutoff=0.65):
    num_atoms = positions.shape[0]
    
    with Pool(4) as p:
        S6 = list(tqdm(p.starmap(process_index, [(i, positions, cutoff) for i in range(num_atoms)]), total=num_atoms, position=0)) 

    return np.array(S6)



if __name__ == '__main__':
    print("Calculating hexagonal order parameters")
    print(average_positions_2d.size)
    print(average_positions_2d.shape)
    S6_last_half_frames = calculate_S6_multiprocessing(average_positions_2d)

    print("Final Frame")
    S6_final_frame = calculate_S6_multiprocessing(final_positions_2d) # for plot
 
    def largest(S6_last_half_frames):
        return reduce(max, S6_last_half_frames)
    
    Ans = largest(S6_last_half_frames)
    print("Largest in given array:", Ans)


    # DEBUG: Check S6 values
    print(f"First few S6 values: {S6_last_half_frames[:10]}")
    print(f"Last few S6 values: {S6_last_half_frames[-10:]}")

    # outside 0 to 1 (DEBUG)
    out_of_range_S6 = S6_last_half_frames[(S6_last_half_frames < 0) | (S6_last_half_frames > 1)]
    if len(out_of_range_S6) > 0:
        print(f"S6 values out of range: {out_of_range_S6}")


    avg_S6 = np.mean(S6_last_half_frames)
    print(f"Average S6 order parameter: {avg_S6}")
    std_S6 = np.std(S6_last_half_frames)
    print(f"Standard deviation of S6 order parameter: {std_S6}")
    sem_S6 = std_S6/np.sqrt(S6_last_half_frames.size)
    print(f"Standard error of S6 order parameter: {sem_S6}")
    print(S6_last_half_frames.size)

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
        scatter = plt.scatter(positions_2d[:, 0], positions_2d[:, 1], c=S6, cmap='coolwarm_r', s=15, vmin=0.0, vmax=1.0, edgecolors=None) # "_r" flips color bar & c=S6 should scale colors...
        plt.colorbar(scatter, label='$S_6$', extend='both')
        
        plt.xlim(np.min(positions_2d[:, 0]), np.max(positions_2d[:, 0])) # x range
        plt.ylim(np.min(positions_2d[:, 1]), np.max(positions_2d[:, 1])) # y range
        plt.gca().set_aspect('equal', 'box') # equal xy scaling
        plt.title(title)

        plt.savefig(title, dpi=300)
        plt.show()

    top_positions = final_positions_2d[:128]  # change to 122 for slipids
    top_S6 = S6_final_frame[:128]
    plot_hexagonal_order_parameter(top_positions, top_S6, 'OpenFF (300K) Top Monolayer Hexagonal Order Parameter ($S_6$)')

    bottom_positions = final_positions_2d[128:]
    bottom_S6 = S6_final_frame[128:]
    plot_hexagonal_order_parameter(bottom_positions, bottom_S6, 'OpenFF (300K) Bottom Monolayer Hexagonal Order Parameter ($S_6$)')
