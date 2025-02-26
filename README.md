# OpenFF-Bilayer

<img width="573" alt="Screenshot 2025-02-26 at 1 31 51 AM" src="https://github.com/user-attachments/assets/afcdc198-a11c-45e1-b87b-d201ccfba8c2" /> </br>

This repository is dedicated to the validation of OpenFF parameters in lipid molecules and bilayer systems. Bilayer construction and parameterization workflows, pure and compositional bilayer systems for OpenFF, and some relevant analysis processes can be found here. </br>

  - Lipid directories contain equilibration and production MD files in addition to topologies before and after refit.
  - Initial bilayer structures were built using OpenFF Interchange pdb structures in Packmol.
  - OpenFF alkane analysis can be found here [[OpenFFLipid](https://github.com/JHoeflich1/OpenFFLipid)].

#### Lipid types with DOIs containing production simulations used for analysis
  1. POPC
     - [2.2.0](10.5281/zenodo.14714284)
     - [Aux](https://zenodo.org/records/14713793)
     - [Core](https://zenodo.org/records/14713704)
  2. POPE
     - [2.2.0](10.5281/zenodo.14714284)
  3. POPS
     - [2.2.0](10.5281/zenodo.14714284)
  4. PSM (SM16)
     - [2.2.0](10.5281/zenodo.14714284)
     - [Aux](https://zenodo.org/records/14713992)
     - [Core](https://zenodo.org/records/14713956)
  5. DMPC
     - [2.2.0](10.5281/zenodo.14714284)
     - [Core](https://zenodo.org/records/14713912)
  6. POPC/15% Cholesterol
     - [2.2.0](10.5281/zenodo.14714284)
     - [Aux](https://zenodo.org/records/14713797)
     - [Core](https://zenodo.org/records/14713823)
  7. POPC/50% Cholesterol
     - [2.2.0](10.5281/zenodo.14714284)
     - [Aux](https://zenodo.org/records/14713870)
     - [Core](https://zenodo.org/records/14713850)

## To-dos
* Add cutoff and water tests.
* Add control descriptions.
* Add Lipid21 systems and CHARMM36 references.
* Zenodo DOIs will contain additional info on OpenFF-specific mods (e.g. cutoff) since that is where it is relevant.
* Add mapping files + NMR ref
* Add a pic
