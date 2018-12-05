#!/home/groups/aboettig/Software/anaconda3/envs/openmm-env/bin/python

import numpy as np
import re
import sys

def read_in_data(filename_list):
        data = list()
        fp_list = open(filename_list, "r")
        searchPattern = re.compile('\s*<Position x.+')

        datasetCount = 0
        for line in fp_list:
                currfile = line.rstrip("\n")

                fp = open(currfile, "r")
                data.append(list())
                for line in fp:
                        if searchPattern.search(line):
                                #print("Match at line [" + line + "]\n")
                                valPattern = re.match(r'.*x="(.+)"\sy="(.*)"\sz="(.*)".*', line)
                                #print("match is [" + valPattern.group(1) + "] [" + valPattern.group(2) + "] [" +valPattern.group(3) +  "]\n")
                                x = float(valPattern.group(1))
                                y = float(valPattern.group(2))
                                z = float(valPattern.group(3))
                                data[datasetCount].append([x,y,z])

                fp.close()
                datasetCount = datasetCount + 1


        fp_list.close()
        return data

def calculate_distance(pt1, pt2):
	return np.sqrt((pt1['x'] - pt2['x'])^2 + (pt1['y'] - pt2['y'])^2 + (pt1['z'] - pt2['z'])^2)


def calculate_radius_of_gyration(xyz):
	rgs = list()	

	for timestep in range(len(xyz)):
		currDat = xyz[timestep]
		normDat = currDat - np.mean(currDat, axis=0)[None,:]
		rgs.append(np.sqrt(np.sum(np.var(np.array(normDat),0))))

	return rgs	

def main():

	if len(sys.argv) < 2:
		print("Usage: calculate_polymer_statistics.py <list_of_state_files.txt> <path/to/outfile/outfile_base_name>\n", file=sys.stderr)
		sys.exit(1)

	xyzData = read_in_data(sys.argv[1])
	outfilename = sys.argv[2]

	rgs = calculate_radius_of_gyration(xyzData)

	fp = open(outfilename + "_RGs.txt", "w")

	for i in range(len(rgs)):
		fp.write(str(rgs[i]) + "\n")

	fp.close()

main()
