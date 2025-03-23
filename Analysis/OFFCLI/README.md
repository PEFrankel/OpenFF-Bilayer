# OpenFF Lipid Bilayer Analysis Tool

A command-line interface for analyzing lipid bilayer molecular dynamics simulations. This tool provides a unified interface for various types of membrane analyses including dihedrals, angles, area per lipid, radial distribution functions, and isomerization behaviors.

+Lateral diffusion measurements, basic MSD plot (with trajectory modification), and torsion potential energy verifications using gmx rerun TBA

## Overview

This collection of Python scripts analyzes molecular dynamics simulations of lipid bilayers. It's designed to work with GROMACS trajectories and provides analyses that are commonly needed in membrane biophysics studies. Most built-in atom naming functionality is intended for use with OpenFF simulations. Isomeric SMILES strings used for reference in the scripts (typically POPC) can be found in the `Simulations` directory in this repository.

## Packages

- Python 3.6+
- MDAnalysis
- NumPy
- Pandas
- Matplotlib
- GROMACS

## Installation

1. Clone the repository:

```bash
git clone https://github.com/yourusername/OpenFF-Bilayer.git
cd Analysis/OFFCLI
```

2. Install required Python packages using your preferred choice of Conda/Mamba, pip, etc.

## Usage

The main script `OFFCLI.py` provides an interface to several analysis scripts:

```bash
python OFFCLI.py --analysis TYPE [OPTIONS]
```

Where `TYPE` can be one of:
- `dihedral`: Analyze dihedral angles
- `angle`: Analyze angles
- `apl`: Calculate area per lipid
- `rdf`: Calculate radial distribution functions
- `isomerization`: Analyze sn-1 gauche isomerization behavior

### Common Options

- Add the following input files to the `/data` directory for better organization:
  - `--xtc`: Path to XTC trajectory file
  - `--tpr`: Path to TPR topology file
  - `--gro`: Path to GRO structure file
- `--output-dir`: Directory for output files (default: "output")

## Output Structure

The results are saved in the specified output directory (default: "output") with the following structure:

```
output/
├── xvg_files/          # Contains XVG formatted data for plots
├── index_files/        # Contains GROMACS index files
├── apl_results*        # Area per lipid results
└── rdf_analysis*       # RDF analysis results
```

## Documentation

Each script has detailed help that can be accessed using the `--help` flag:

```bash
python OFFCLI.py --help
python OFFCLI.py --analysis dihedral --help
```

## Example Commands

### Dihedral Analysis

```bash
python OFFCLI.py --analysis dihedral \
  --xtc path/to/trajectory.xtc \
  --nlip 128 \
  --offset 134 \
  --torsion_set sn1_isomerization
```

This analyzes the dihedral angles related to sn-1 chain isomerization for 128 lipids, each containing 134 atoms.

### Radial Distribution Function Analysis

```bash
python OFFCLI.py --analysis rdf \
  --xtc path/to/trajectory.xtc \
  --tpr path/to/topology.tpr
```

This calculates the radial distribution function between phosphate and amine groups across lipids. This avoids conflicts on the same lipid to prevent overestimation of head group interactions.

### Area Per Lipid Analysis

```bash
python OFFCLI.py --analysis apl \
  --xtc path/to/trajectory.xtc \
  --gro path/to/structure.gro \
  --resname POPC \
  --max-time 500000
```

This example calculates the area per lipid for POPC lipids up to 500 ns.

### Isomerization Analysis

```bash
python OFFCLI.py --analysis isomerization \
  --title "Force Field Comparison"
```

This script compares isomerization data from different force fields, assuming the necessary XVG files are already generated.