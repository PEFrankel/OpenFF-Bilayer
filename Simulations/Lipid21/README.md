## GROMACS-Converted Lipid21 systems
- https://zenodo.org/records/14776136

## Workflow for generating Lipid21 systems for GROMACS

1. Download a version of AMBER: https://ambermd.org/GetAmber.php
2. Download acpype for conversion: https://github.com/alanwilter/acpype
3. https://github.com/callumjd/lipid21 (this repository has some build instructions too)
    - pdb4amber if using a structure built on your own and it's not already initialized
4. $ tleap -s -f leaprc.lipid21 
    - (in the relevant directory with the PDB mentioned below)
5. Leap script using bilayer pdb
6. $ acpype -p POPC_128.prmtop -x POPC_128.inpcrd
    - Converts amber to gmx
