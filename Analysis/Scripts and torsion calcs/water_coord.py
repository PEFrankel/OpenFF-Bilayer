import mdtraj as md
import numpy as np

traj = md.load('md.xtc', top='md.gro')
headgroup_indices = traj.topology.select("resname POPC and name P1x O1x O2x O3x O4x N1x C3x C4x C5x C6x C7x H4x H5x H6x H7x H8x H9x H10x H11x H12x H13x H14x H15x H16x")
water_indices = traj.topology.select("resname TIP3P and name O1w H1w H2w")
# Slipids
# headgroup_indices = traj.topology.select("resname POPC and name P O11 O12 O13 O14 N C11 C12 C13 C14 C15 H15A H15B H15C H14A H14B H14C H13A H13B H13C H12A H12B H11A H11B")
# water_indices = traj.topology.select("resname SOL and name OH2 H1 H2")
# MacRog
# headgroup_indices = traj.topology.select("resname POPC and name P1 O1 O2 O3 O4 N1 C1 C2 C3 C4 C5 H1 H2 H3 H4 H5 H6 H7 H8 H9 H10 H11 H12 H13")
# water_indices = traj.topology.select("resname SOL and name OW HW1 HW2")

# DEBUG
print(f"Number of headgroup atoms selected: {len(headgroup_indices)}")
print(f"Number of water oxygen atoms selected: {len(water_indices)}")

if len(headgroup_indices) == 0 or len(water_indices) == 0:
    raise ValueError("One of the selections is empty. Check the atom selection criteria.")

cutoff_distance = 0.35
last_frame = traj[-1000]

atom_pairs = np.array([[h, w] for h in headgroup_indices for w in water_indices], dtype=np.int32)
distances = md.compute_distances(last_frame, atom_pairs)
distances = distances.reshape(len(headgroup_indices), len(water_indices))
num_coordinated_waters = np.sum(distances < cutoff_distance, axis=1)

avg_coordination = np.mean(num_coordinated_waters)
std_coordination = np.std(num_coordinated_waters)
print(f"Average Coordination Number in the last frame: {avg_coordination:.3f}")
print(f"Standard Deviation of Coordination Number: {std_coordination:.3f}")
