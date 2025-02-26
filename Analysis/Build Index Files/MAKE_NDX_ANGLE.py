name = "O2_P1_O4"
fname = name + ".ndx"
label = "[ " + name + " ]"
nlip = 122
offset = 134
atoms = [22,20,23] # "O2P_P_O2P"
pfile = open(fname,"w")
print(label,file=pfile)
for n in range(nlip):
    print(f"{atoms[0]+offset*n:4d} {atoms[1]+offset*n:4d} {atoms[2]+offset*n:4d}",file=pfile)
pfile.close()
