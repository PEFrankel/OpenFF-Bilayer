# OpenFF_Bilayer

OpenFF pure and compositional bilayers.
  Lipid types with DOIs containing trajectories and topologies used for analysis:
  1. POPC (https://zenodo.org/records/10265709)

    - TIP3P waters
  2. POPE (https://zenodo.org/records/10265965)

    - 310K
    - TIP3P waters
  3. POPS (https://zenodo.org/records/10266049)

    - Na+ counterions
    - SPC/E waters
  4. POPC/Cholesterol (https://zenodo.org/records/10344232)

    - 20% chol
    - SPC/E waters
  5. POPC/Cholesterol (https://zenodo.org/records/10366462)

    - 15% chol
    - SPC/E waters

  - Trajectories contain the 200 ns MD production simulations with data saved every 10 ps.
  - Lipid directories contain equilibration and production MD files.
  - Directories also contain parametrization of lipid molecules using the OpenFF Toolkit and Interchange conversion method. Residue names were changed to assemble the final topology.
  - The Methods directory contains commands used for running the simulations in GROMACS as well as a script for building the 15% cholesterol bilayer.

  - The PDB for POPS was sourced from: http://www.fos.su.se/~sasha/SLipids/Downloads.html. This was then converted to an SDF for use in the OpenFF Toolkit using OpenBabel.
  - All initial bilayer gro files except for 15% cholesterol were adopted from simulations submitted to the NMRLipids Databank and reformatted for naming continuity.
  - Bilayers were equilibrated for 100 ns.
  - Systems were determined to be equilibrated by the plateau of instantaneous temperature, pressure, density, and potential energy.
  - Further analysis of equilibration using PCA is currently being worked on.

A control using MacRog forcefields and OpenFF methods was run and compared to its existing simulation analysis in the NMRLipids Databank. MacRog was chosen for its similar nomenclature and atom ordering for ease of topology creation.

MacRog POPC DOI: https://zenodo.org/record/3741793

MacRog POPE DOI: https://zenodo.org/record/3725670
