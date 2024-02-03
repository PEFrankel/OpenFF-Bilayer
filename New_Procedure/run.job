#!/bin/bash
#SBATCH -N 1 --ntasks-per-node=64
#SBATCH -t 24:00:00
#SBATCH -p amilan
#SBATCH -J 'POPE'
#SBATCH -o '%x.out'
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pafr7911@colorado.edu

# load in necessary modules (these are the ones for Mando's Gromacs installation)
ml gcc/11.2.0
ml openmpi/4.1.1

# source Gromacs
source /projects/dora1300/pkgs/gromacs-2022-cpu-mpi/bin/GMXRC

# Energy minimizaton
#mpirun -np 1 gmx_mpi grompp -p topol.top -f min.mdp -c bilayer.gro -o min.tpr 
#mpirun -np 64 gmx_mpi mdrun -deffnm min

# NVT equilibration (40 ns)
#mpirun -np 1 gmx_mpi grompp -p topol.top -f nvt.mdp -c min.gro -o nvt.tpr
mpirun -np 64 gmx_mpi mdrun -deffnm nvt

# NPT equilibration (60 ns)
mpirun -np 1 gmx_mpi grompp -p topol.top -f npt.mdp -c nvt.gro -o npt.tpr
mpirun -np 64 gmx_mpi mdrun -deffnm npt
#mpirun -np 64 gmx_mpi mdrun -deffnm npt -cpi npt.cpt

# Production MD
mpirun -np 1 gmx_mpi grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
mpirun -np 64 gmx_mpi mdrun -deffnm md
#mpirun -np 64 gmx_mpi mdrun -deffnm md -cpi md.cpt
