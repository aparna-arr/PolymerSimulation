# MAKE SURE THIS IS RUN UNDER SYSTEM'S PYTHON3, NOT ANACONDA
# or else you can't import anything

from openmmlib import polymerutils
import sys
import os
import subprocess

if len(sys.argv) < 5:
	print("usage: generate_polymer_length_dist.py <start of range: 2^x> <end of range: 2^y> </path/to/output/dir> <out filename prefix>", file=sys.stderr)
	sys.exit(1)

startExp = int(sys.argv[1])
endExp = int(sys.argv[2])
outDir = sys.argv[3]
outPrefix = sys.argv[4]

if not os.path.isdir(outDir):
	subprocess.run(["mkdir", "-p", outDir])

fh_index = open(outDir + "/" + outPrefix + "_contents.txt", "w")

for k in range(startExp,endExp+1):
	polymer_len = 2**k

	print("startExp: [", startExp, "] endExp: [", endExp, "] k [", k, "] polymer_len [",polymer_len,"]")
	
	polymer = polymerutils.create_spiral(r1=10, r2=13, N=polymer_len)
	polymer = polymer.transpose()

	filename = outDir + "/" + outPrefix + "_" + str(polymer_len) + ".xyz";
	fh = open(filename, "w");

	fh.write(str(len(polymer)))
	fh.write("\n\n")

	for i in range(len(polymer)):
		fh.write("A\t" + str(polymer[i,0]) + "\t" + str(polymer[i,1]) + "\t" + str(polymer[i,2]) + "\n")

	fh.close()

	fh_index.write(filename + "\n")

fh_index.close()
