import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi
from multiprocessing import Pool, cpu_count
from tqdm import tqdm
import time

# Load trajectory and topology
trajectory = md.load('OpenFF_half_mol.xtc', top='OpenFF.gro')

# Define acyl chain patterns
acyl_chain_patterns = [
    ['C8x', 'C9x', 'C10x', 'C11x', 'C12x', 'C13x', 'C14x', 'C15x', 'C16x', 'C17x', 'C18x', 'C19x', 'C20x', 'C21x', 'C22x', 'C23x', 'C24x', 'C25x'],
    ['C27x', 'C28x', 'C29x', 'C30x', 'C31x', 'C32x', 'C33x', 'C34x', 'C35x', 'C36x', 'C37x', 'C38x', 'C39x', 'C40x', 'C41x', 'C42x']
]

def separate_acyl_chains(topology, patterns):
    chains = []
    for residue in topology.residues:
        for pattern in patterns:
            chain = [atom.index for atom in residue.atoms if atom.name in pattern]
            if len(chain) == len(pattern):
                chains.append(chain)
    return chains

acyl_chains = separate_acyl_chains(trajectory.topology, acyl_chain_patterns)

def compute_frame_averages(positions, chains):
    print("Computing frame averages...")
    frame_averages = []
    for frame_idx, frame_positions in enumerate(positions):
        frame_avg_positions = []
        for chain in chains:
            chain_positions = frame_positions[chain, :2]  # Get x and y positions for the current chain
            avg_position = np.mean(chain_positions, axis=0)
            frame_avg_positions.append(avg_position)
        frame_averages.append(np.array(frame_avg_positions))
    print("Frame averages computed.")
    return np.array(frame_averages)

def calculate_S6_for_chain(i, positions, cutoff=0.65):
    delta = positions - positions[i]
    distances = np.sqrt(np.sum(delta**2, axis=1))
    within_cutoff = (distances < cutoff)
    within_cutoff[i] = False
    neighbor_indices = np.where(within_cutoff)[0]

    if len(neighbor_indices) == 0:
        return 0.0

    angles = []
    for neighbor in neighbor_indices:
        delta = positions[neighbor] - positions[i]
        angle = np.arctan2(delta[1], delta[0])
        angles.append(angle)

    S6 = (1/6) * np.sqrt(
        (np.sum(np.cos(6 * np.array(angles)))**2 +
         np.sum(np.sin(6 * np.array(angles)))**2)
    )
    return S6

def process_frame(frame_idx, positions, chains, cutoff=0.65):
    print(f"Processing frame {frame_idx}...")
    frame_positions = positions[frame_idx]
    S6_results = []
    for chain in chains:
        chain_positions = frame_positions[chain, :2]  # Get x and y positions for the current chain
        for i in range(len(chain_positions)):
            S6 = calculate_S6_for_chain(i, chain_positions, cutoff)
            S6_results.append(S6)
    print(f"Frame {frame_idx} processed.")
    return S6_results

def calculate_S6_across_frames(positions, chains, cutoff=0.65):
    num_frames = positions.shape[0]
    print(f"Number of frames in trajectory: {num_frames}")
    all_S6_results = []

    # Use multiprocessing to process frames in parallel
    with Pool(cpu_count()) as pool:
        results = []
        for frame_idx in range(num_frames):
            results.append(pool.apply_async(process_frame, (frame_idx, positions, chains, cutoff)))
        
        for i, result in enumerate(results):
            frame_results = result.get()
            all_S6_results.append(frame_results)
            print(f"Frame {i} results collected.")
    
    return np.array(all_S6_results)

def compute_statistics(S6_results):
    print("Computing aggregate statistics...")
    mean_S6 = np.mean(S6_results, axis=0)
    std_S6 = np.std(S6_results, axis=0)
    sem_S6 = std_S6 / np.sqrt(S6_results.shape[0])
    print("Aggregate statistics computed.")
    return mean_S6, std_S6, sem_S6

if __name__ == '__main__':
    positions = trajectory.xyz  # Get all frames' positions
    print(f"Loaded {positions.shape[0]} frames from trajectory.")
    
    frame_averages = compute_frame_averages(positions, acyl_chains)

    # Split data for multiprocessing
    print("Starting S6 calculations...")
    S6_results = calculate_S6_across_frames(positions, acyl_chains)
    print("S6 calculations completed.")

    mean_S6_frames, std_S6_frames, sem_S6_frames = compute_statistics(S6_results)
    mean_S6 = np.mean(mean_S6_frames)
    std_S6 = np.mean(std_S6_frames)
    sem_S6 = np.mean(sem_S6_frames)

    print(f"Mean S6: {mean_S6}")
    print(f"STD S6: {std_S6}")
    print(f"SEM S6: {sem_S6}")

    # Plot the final frame
    final_frame_positions = positions[-1]
    final_frame_avg_positions = compute_frame_averages(final_frame_positions[np.newaxis, :], acyl_chains)[0]

    # plt.figure(figsize=(10, 10))
    # plt.scatter(final_frame_avg_positions[:, 0], final_frame_avg_positions[:, 1], s=10)
    # plt.title('Final Frame Acyl Chain Positions')
    # plt.xlabel('X Position')
    # plt.ylabel('Y Position')
    # plt.show()
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

    top_positions = final_frame_avg_positions[:128]  # change to 122 for slipids
    top_S6 = final_frame_avg_positions[:128]
    plot_hexagonal_order_parameter(top_positions, top_S6, 'OpenFF (300K) Top Monolayer Hexagonal Order Parameter ($S_6$)')

    bottom_positions = final_frame_avg_positions[128:]
    bottom_S6 = final_frame_avg_positions[128:]
    plot_hexagonal_order_parameter(bottom_positions, bottom_S6, 'OpenFF (300K) Bottom Monolayer Hexagonal Order Parameter ($S_6$)')
