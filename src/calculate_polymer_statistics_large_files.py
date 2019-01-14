#!/home/groups/aboettig/Software/anaconda3/envs/openmm-env/bin/python

import numpy as np
import re
import sys

def calculate_rgs(filename_list, polymer_len, domains):
	print("In function calculate_rgs()", flush=True)

	print("Opening state file list", flush=True)
	fp_list = open(filename_list, "r")
	searchPattern = re.compile('\s*<Position x.+')

	state_file_list = list()	

	for line in fp_list:
		currfile = line.rstrip("\n")
		state_file_list.append(currfile)
	fp_list.close()

	print("Done reading state file list", flush=True)

	numTimepoints = len(state_file_list)
	rgs = list()	

	print("Begin loop over state files", flush=True)
	for currfile in state_file_list:
		print("Opening file[" + currfile + "]", flush=True)

		fp = open(currfile, "r")
		data = np.zeros((polymer_len,3))
		lineCount = 0

		for line in fp:
			if searchPattern.search(line):
				#print("Match at line [" + line + "]\n")
				valPattern = re.match(r'.*x="(.+)"\sy="(.*)"\sz="(.*)".*', line)
				#print("match is [" + valPattern.group(1) + "] [" + valPattern.group(2) + "] [" +valPattern.group(3) +  "]\n")
				x = float(valPattern.group(1))
				y = float(valPattern.group(2))
				z = float(valPattern.group(3))

				data[lineCount,0] = x
				data[lineCount,1] = y
				data[lineCount,2] = z		
				#print("lineCount is [" + str(lineCount) + "]", flush=True)	
				lineCount = lineCount + 1	

		fp.close()
		
		print("Closing file [" + currfile + "]", flush=True)		

		#print("Normalizing data", flush=True)

		normDat = data - np.mean(data, axis=0)[None,:]
		
		#print("Done normalizing data", flush=True)		

		if len(domains) == 0:
			print("Calculating RGS: no domains", flush=True)
			rgs.append(str(np.sqrt(np.sum(np.var(np.array(normDat),0)))))
			print("Done calculating RGS: no domains", flush=True)
		else:
			print("Calculating RGS: domains (full polymer)", flush=True)
			domainRgsStr = str(np.sqrt(np.sum(np.var(np.array(normDat),0))))
			print("Done calculating RGS: domains (full polymer)", flush=True)
			for dStart,dEnd in domains:
				#print("Calculating rgs for domain [" + str(dStart) + "," + str(dEnd) + "]", flush=True)
				#print("Extracting data", flush=True)
				domainDat = data[dStart:dEnd]
				#print("Data shape: [" + str(domainDat.shape) +"]", flush=True)
				#print("Data is [" + str(domainDat) + "]",flush=True)
				#print("Normalizing data", flush=True)
				domainNormDat = domainDat - np.mean(domainDat, axis=0)[None,:]
				#print("Calculating RGS", flush=True)
				domainRgs = str(np.sqrt(np.sum(np.var(np.array(domainNormDat),0))))
				domainRgsStr = domainRgsStr + '\t' + domainRgs	
				#print("Done with this domain", flush=True)
			
			rgs.append(domainRgsStr)

	print("Done with state file loop", flush=True)
	return rgs

def main():
	print("In Main", flush=True)

	if len(sys.argv) < 3:
		print("Usage: calculate_polymer_statistics_large_files.py <list_of_state_files.txt> <polymer length> <path/to/outfile/outfile_base_name> <optional: comma-sep specific domains to calculate rgs within (0-based) like 1000,2000 OR a single integer representing step size to calculate RGs in windows across polymer>\n", file=sys.stderr)
		sys.exit(1)

	print("Getting arguments", flush=True)

	statefileList = sys.argv[1]
	polymer_len = int(sys.argv[2])
	outfilename = sys.argv[3]
	domains = list()
	header = "FullPolymer"

	if len(sys.argv) > 4:
		if ',' in sys.argv[4]:
			for i in range(4,len(sys.argv)):
				start,end = sys.argv[i].split(',')
				domains.append((int(start),int(end)))
				header = header + '\t' 'domain_' + start + "_" + end
		else:
			stepsize = int(sys.argv[4])
			for i in range(0,polymer_len,stepsize):
				end = i + stepsize
				if end > polymer_len:
					end = polymer_len

				domains.append((i,end))
				header = header + '\t' + 'domain_' + str(i) + "_" + str(end)

	print("Calculating RGS", flush=True)
	rgs = calculate_rgs(statefileList, polymer_len, domains)
	print("Done calculating RGS", flush=True)

	fp = open(outfilename + "_RGs.txt", "w")
	fp.write(header + "\n")

	for i in range(len(rgs)):
		fp.write(rgs[i] + "\n")

	fp.close()

print("Starting script, calling Main()", flush=True)

main()
