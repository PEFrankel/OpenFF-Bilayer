# OpenFF_Bilayer

OpenFF pure POPC and POPE bilayers. 128 lipid systems hydrated with 5000 TIP3P waters. Standard conditions: 1 Bar, 300 K, and no ions. Methods directories contain equilibration files as well as relevant trajectories and topologies used for analysis. Trajectories contain the 200 ns MD production simulations with data saved every 10 ps. Parameter directories contain parametrization of POPC and POPE molecules using the OpenFF Toolkit and Interchange conversion method. This simulation will be submitted to the NMRLipid Databank respository and references will be added to Table S6 of the manuscript:

https://github.com/NMRLipids/DataBankManuscriptText/blob/main/SI.pdf

A control using MacRog forcefields and OpenFF methods was run and compared to its existing simulation analysis in the NMRLipids Databank. MacRog was chosen for its similar nomenclature and atom ordering for ease of topology creation.

MacRog POPC DOI: https://zenodo.org/record/3741793
