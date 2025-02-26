# POPC/50% CHOL

Lipids: POPC/CHOL<br/>
Lipid count: 64 POPC (64 per leaflet) / 64 CHOL (64 per leaflet)<br/>
Solvent: TIP3P<br/>
Solvent count: 5120<br/>
Ions: no<br/>
Ion Count: n/a<br/>
Temperature: 300K<br/>
Length of equilibration: 100ns (40 NVT - 60 NPT)<br/>
Length of production: 500ns<br/>

- HMR w/ 4fs timestep
- 0.9 rvdw/rcouloumb cutoff
- 0.8 rvdw-switch

SMILES strings: <br/>
POPC: `[C@](COP(=O)([O-])OCC[N+](C)(C)C)([H])(OC(CCCCCCC/C=C\CCCCCCCC)=O)COC(CCCCCCCCCCCCCCC)=O`<br/>
CHOL: `[C@]12(CC=C3C[C@@H](O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@@]([H])([C@@](C)([H])CCCC(C)C)CC[C@@]21[H])[H]`