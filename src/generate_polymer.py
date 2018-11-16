# MAKE SURE THIS IS RUN UNDER SYSTEM'S PYTHON3, NOT ANACONDA
# or else you can't import anything

from openmmlib import polymerutils

polymer = polymerutils.create_spiral(r1=10, r2=13, N=5000)
polymer = polymer.transpose()

filename = "testpolymer_spiral.xyz";
fh = open(filename, "w");

fh.write(str(len(polymer)))
fh.write("\n\n")

for i in range(len(polymer)):
	fh.write("A\t" + str(polymer[i,0]) + "\t" + str(polymer[i,1]) + "\t" + str(polymer[i,2]) + "\n")

fh.close()
