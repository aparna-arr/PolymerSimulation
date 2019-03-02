#!/share/software/user/open/python/3.6.1/bin/python3
#module load python/3.6.1
# MAKE SURE THIS IS RUN UNDER SYSTEM'S PYTHON3, NOT ANACONDA
# or else you can't import anything

from openmmlib import polymerutils
import sys
import os
import subprocess
import math

def read_input_file(infile):
	params = dict()
	params['PARAM_FILE'] = infile
	params['DOMAIN_LIST'] = list()

	fp = open(infile, 'r')
	
	blockFlag = False
	currBlockIdx = 0
	for line in fp:
		if line.startswith('#') or not line.strip():
			continue

		paramAr = line.rstrip().split('\t')
		#print("line is [" + line + "]")
		paramKey = paramAr[0]
	
		#print("paramKey is [" + paramKey + "]")	
		if paramKey.startswith('@') or blockFlag == True:
			# in a DOMAIN block
			if paramKey == "@END":
				blockFlag = False
				currBlockIdx = currBlockIdx + 1
			elif blockFlag == True:
				#paramVal = paramAr[1].split('#',1)[0].strip()
				paramVal = paramAr[1]
				if paramKey in ["DOMAIN_MARKER"]:
					params['DOMAIN_LIST'][currBlockIdx][paramKey] = paramVal
				elif paramKey in ["MARKER_PERC"]:
					params['DOMAIN_LIST'][currBlockIdx][paramKey] = float(paramVal)
				elif paramKey in ['DOMAIN_START', 'DOMAIN_END']:
					params['DOMAIN_LIST'][currBlockIdx][paramKey] = int(paramVal)
			elif paramKey == "@DOMAIN":
				params['DOMAIN_LIST'].append(dict())
				blockFlag = True
		else:
			paramVal = paramAr[1].split('#',1)[0].strip()
			if paramKey in ["POLYMER_LEN"]:
				params[paramKey] = int(paramVal)	
			elif paramKey in ["DEFAULT_MARKER", "SAVE_PATH", "SAVE_FILENAME"]:
				params[paramKey] = paramVal
			

	posToIdx = dict()

	for i in range(len(params['DOMAIN_LIST'])):
		start = params['DOMAIN_LIST'][i]['DOMAIN_START']
		posToIdx[start] = i 

	return params,posToIdx	


def main():
	if len(sys.argv) < 2:
		print("usage: generate_polymer_domain.py <domain specification input file>", file=sys.stderr)
		sys.exit(1)

	params,posToIdx = read_input_file(sys.argv[1])

	if not os.path.isdir(params['SAVE_PATH']):
		subprocess.run(["mkdir", "-p", params['SAVE_PATH']])

	#fh_index = open(params['SAVE_PATH'] + "/" + params['SAVE_FILENAME'] + "_contents.txt", "w")

	polymer_len = params['POLYMER_LEN']

	polymer = polymerutils.create_spiral(r1=10, r2=13, N=polymer_len)
	polymer = polymer.transpose()

	filename = params['SAVE_PATH'] + "/" + params['SAVE_FILENAME']
	con = open(params['SAVE_PATH'] + "/" + "contents.txt", "w")
	con.write(filename)
	con.close()

	fh = open(filename, "w");

	fh.write(str(len(polymer)))
	fh.write("\n\n")

	blockFlag = False
	blockSpacer = 1
	blockMarker = params['DEFAULT_MARKER']
	currDomainStart = 0
	for i in range(len(polymer)):
		letter = params['DEFAULT_MARKER']

		if blockFlag == True:
			if i >= params['DOMAIN_LIST'][posToIdx[currDomainStart]]['DOMAIN_END']:
				blockFlag = False
				blockSpacer = 1
				blockMarker = params['DEFAULT_MARKER']
				currDomainStart = len(posToIdx)
			else:
				# turn to block letter depending on percent, use some modulo trick
				if i % blockSpacer == 0:
					letter = blockMarker
		elif i in posToIdx:
			blockFlag = True
			currDomainStart = i
			# turn to block letter depending on percent, use some modulo trick
			perc = params['DOMAIN_LIST'][posToIdx[i]]['MARKER_PERC']
			blockMarker = params['DOMAIN_LIST'][posToIdx[i]]['DOMAIN_MARKER']
			blockSpacer = round(1 / perc,0)

			#print("blockMarker [" + blockMarker + "]")
		
			#print("i is [" + str(i) + "]")
			#print("blockSpacer [" + str(blockSpacer) + "]")
			#print("i % blockSpacer is [" + str(i % blockSpacer) + "]")

			if i % blockSpacer == 0:
				letter = blockMarker

		fh.write(letter + "\t" + str(polymer[i,0]) + "\t" + str(polymer[i,1]) + "\t" + str(polymer[i,2]) + "\n")

	fh.close()

	#fh_index.write(filename + "\n")

	#fh_index.close()

main()
