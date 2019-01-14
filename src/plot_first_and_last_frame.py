#!/home/groups/aboettig/Software/anaconda3/envs/openmm-env/bin/python

import numpy as np
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.animation as animation
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import re
import math
import matplotlib.cm as cm


def read_in_data(filename_list):
	data = list()

	max_x = 0
	min_x = 0
	max_y = 0
	min_y = 0
	max_z = 0
	min_z = 0

	fp_list = open(filename_list, "r")
	searchPattern = re.compile('\s*<Position x.+')
	
	datasetCount = 0
	lineCount = 0
	fileList = list()
	currfile = ""
	for line in fp_list:
		currfile = line.rstrip("\n")
		
		if lineCount == 0:
			fileList.append(currfile)

		lineCount = lineCount + 1

	if lineCount > len(fileList) + 1:
		fileList.append(currfile)
	
	for f in fileList:
		fp = open(f, "r")
		data.append(dict())
		data[datasetCount] = {
			'x' : list(),
			'y' : list(),
			'z' : list()
		}
		for line in fp:
			if searchPattern.search(line):
				#print("Match at line [" + line + "]\n")
				valPattern = re.match(r'.*x="(.+)"\sy="(.*)"\sz="(.*)".*', line)
				#print("match is [" + valPattern.group(1) + "] [" + valPattern.group(2) + "] [" +valPattern.group(3) +  "]\n")
				x = float(valPattern.group(1))
				y = float(valPattern.group(2))
				z = float(valPattern.group(3))

				if x < min_x:
					min_x = x
				elif x > max_x:
					max_x = x

				if y < min_y:
					min_y = y
				elif y > max_y:
					max_y = y
			
				if z < min_z:
					min_z = z
				elif z > max_z:
					max_z = z

				data[datasetCount]['x'].append(x)
				data[datasetCount]['y'].append(y)
				data[datasetCount]['z'].append(z)

		fp.close()
		#data[datasetCount] = np.array(data[datasetCount])

		#print("Shape for dataset [" + str(datasetCount) + "] is [" + str(data[datasetCount].shape) + "]\n")


		datasetCount = datasetCount + 1


	fp_list.close()

	#print("minx: " + str(min_x) + " min_y: " + str(min_y) + " min_z: " + str(min_z))
	#print("maxx: " + str(max_x) + " maxy: " + str(max_y) + " maxz: " + str(max_z))

	minmaxes = {
		'min_x' : math.floor(min_x),	
		'min_y' : math.floor(min_y),	
		'min_z' : math.floor(min_z),	
		'max_x' : math.ceil(max_x),
		'max_y' : math.ceil(max_y),
		'max_z' : math.ceil(max_z)
	
	}
	return data, minmaxes

def main():
	
	if len(sys.argv) < 2:
		print("Usage: plot_first_and_last_frame.py <list_of_state_files.txt> <path/to/outfile/outfile_base_name>\n", file=sys.stderr)
		sys.exit(1)

	xyzData, minmax = read_in_data(sys.argv[1])
	outfilename = sys.argv[2]

	numPoints = len(xyzData[0]['x'])

	cmap = cm.rainbow(np.linspace(0.0, 1.0, numPoints))

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	first = ax.scatter([], [], [], color=cmap[range(numPoints)])
	first._offsets3d = (xyzData[0]['x'], xyzData[0]['y'], xyzData[0]['z'])
	title = ax.set_title('Polymer Frame First')
	ax.set_xlim(minmax['min_x'], minmax['max_x'])
	ax.set_ylim(minmax['min_y'], minmax['max_y'])
	ax.set_zlim(minmax['min_z'], minmax['max_z'])
	
	fig.savefig(outfilename + '_first.png', bbox_inches='tight')

	fig2 = plt.figure()
	ax2 = fig2.add_subplot(111, projection='3d')
	last = ax2.scatter([], [], [], color=cmap[range(numPoints)])
	last._offsets3d = (xyzData[1]['x'], xyzData[1]['y'], xyzData[1]['z'])
	title2 = ax2.set_title('Polymer Frame Last')
	
	ax2.set_xlim(minmax['min_x'], minmax['max_x'])
	ax2.set_ylim(minmax['min_y'], minmax['max_y'])
	ax2.set_zlim(minmax['min_z'], minmax['max_z'])

	fig2.savefig(outfilename + '_last.png', bbox_inches='tight')
main()
