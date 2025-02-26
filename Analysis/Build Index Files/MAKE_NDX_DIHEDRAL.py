name = "2-5_HC"
fname = name + ".ndx"
label = "[ " + name + " ]"
nlip = 122  #check system lipid count
offset = 134 #atoms per lipid
atoms = [34,33,45,48] # first indices in GRO files (increments using offset)
pfile = open(fname,"w")
print(label,file=pfile)
for n in range(nlip):
    print(f"{atoms[0]+offset*n:4d} {atoms[1]+offset*n:4d} {atoms[2]+offset*n:4d} {atoms[3]+offset*n:4d}",file=pfile)
pfile.close()
