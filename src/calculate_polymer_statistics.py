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


def calculate_radius_of_gyration(xyz, domains):
	rgs = list()	

	for timestep in range(len(xyz)):
		currDat = xyz[timestep]
		normDat = currDat - np.mean(currDat, axis=0)[None,:]
		
		if len(domains) == 0:
			rgs.append(str(np.sqrt(np.sum(np.var(np.array(normDat),0)))))
		else:
			domainRgsStr = str(np.sqrt(np.sum(np.var(np.array(normDat),0))))
			for dStart,dEnd in domains:
				domainDat = currDat[dStart:dEnd]
				domainNormDat = domainDat - np.mean(domainDat, axis=0)[None,:]
				domainRgs = str(np.sqrt(np.sum(np.var(np.array(domainNormDat),0))))
				domainRgsStr = domainRgsStr + '\t' + domainRgs	
			
			rgs.append(domainRgsStr)

	return rgs	

def main():
	print("In Main", flush=True)

	if len(sys.argv) < 3:
		print("Usage: calculate_polymer_statistics.py <list_of_state_files.txt> <path/to/outfile/outfile_base_name> <optional: comma-sep specific domains to calculate rgs within (0-based) like 1000,2000>\n", file=sys.stderr)
		sys.exit(1)

	print("Reading in data", flush=True)
	xyzData = read_in_data(sys.argv[1])
	print("Data read in", flush=True)
	outfilename = sys.argv[2]

	domains = list()
	header = "FullPolymer"
	if len(sys.argv) > 3:
		for i in range(3,len(sys.argv)):
			start,end = sys.argv[i].split(',')
			domains.append((int(start),int(end)))
			header = header + '\t' 'domain_' + start + "_" + end

	print("Calculating RGS", flush=True)
	rgs = calculate_radius_of_gyration(xyzData, domains)
	print("Done calculating RGS", flush=True)

	fp = open(outfilename + "_RGs.txt", "w")
	fp.write(header + "\n")
	for i in range(len(rgs)):
		fp.write(rgs[i] + "\n")

	fp.close()

print("Starting script, calling Main()", flush=True)
main()
