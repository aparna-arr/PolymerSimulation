import numpy as np
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.animation as animation
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import re
import math

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

	for line in fp_list:
		currfile = line.rstrip("\n")

		fp = open(currfile, "r")
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
		print("Usage: plot_animation.py <list_of_state_files.txt> <path/to/outfile/outfile_base_name>\n", file=sys.stderr)
		sys.exit(1)

	xyzData, minmax = read_in_data(sys.argv[1])
	outfilename = sys.argv[2]
	nfr = len(xyzData)
	fps = 10 # Frame per sec
	#print("frames is " + str(nfr))

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	#sct, = ax.plot([],[],[], "bo-", markersize=2)
	#plot(x, y, color='green', linestyle='dashed', marker='o',
	#markerfacecolor='blue', markersize=12).

	sct, = ax.plot([], [], [], color='red', linestyle='solid', marker='.', markerfacecolor='black', markeredgecolor='black', markersize=2)

	def updateAni(ifrm):
		print('ifrm is [', ifrm, ']')
		sct.set_data(xyzData[ifrm]['x'], xyzData[ifrm]['y'])
		sct.set_3d_properties(xyzData[ifrm]['z'])
		return sct
	
	ax.set_xlim(minmax['min_x'], minmax['max_x'])
	ax.set_ylim(minmax['min_y'], minmax['max_y'])
	ax.set_zlim(minmax['min_z'], minmax['max_z'])
	ani = animation.FuncAnimation(fig,updateAni,nfr, interval=50)

	ani.save(outfilename + '.gif', writer='imagemagick', fps = fps)
	
main()
