## Use these steps to produce interpretable input files for the S6 code

1. Get trr and tpr
2. Generate a full gro file from a trr using `-pbc nojump`: **Keep track of how many frames are in this and adjust start time as needed:<br/>
```$ gmx trjconv -f md.trr -s md.tpr -o traj_new.gro -b 100000 -pbc mol (no -n index.ndx?)```<br/>
```$ Select just lipids when prompted```
3. Make a gro file for 1 frame:<br/>
```$ gmx trjconv -f md.trr -s md.tpr -o traj_1f.gro -b 500000 -pbc mol``` **You can specify the last frame in your trajectory for this
4. Generate an index file using the 1 frame gro file with `GENERATE_NDX.py` in Fortran directory
5. Adjust the `bilayer_order_s6.f9` code to recognize your input file names and number of trajectory frames in your .gro
6. Create a mamba environment (or similar) with gfortran installed:<br/>
```$ gfortran bilayer_order_s6.f95 -o order.exe```
7. Run with echo specifying the first three inputs: traj_new.gro (trajectory from step 3), 401 (how many frames are in ‘traj_new.gro’), ndxfile.txt (index file from step 5):<br/>
```$ echo traj_new.gro 401 ndxfile.txt s6.xvg allmaps.xyz hex.xyz | ./order.exe```
