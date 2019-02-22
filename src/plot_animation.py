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
	currLineCount = 0
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

def getColors(domainFile):
	#cmap = cm.rainbow(np.linspace(0.0, 1.0, numPoints))
	colors = list()
	fp = open(domainFile, 'r')
	prev_marker = "INIT"
	prev_count = 0
	total_count = 0
	for line in fp:
		lineAr = line.rstrip().split('\t')
		if len(lineAr) < 4:
			continue

		total_count = total_count + 1
		marker = lineAr[0]
		if marker != prev_marker and prev_marker != "INIT":
			#FIXME hardcoding!
			if prev_marker == "A":
				print("Adding rainbow colors")
				colors.extend(cm.rainbow(np.linspace(0.0, 1.0, prev_count)))
			elif prev_marker == "B":
				print("adding bone colors")
				colors.extend(cm.bone(np.linspace(0.0, 1.0, prev_count)))
			else:
				print("adding copper colors")
				colors.extend(cm.copper(np.linspace(0.0, 1.0, prev_count)))


			prev_marker = marker
			prev_count = 1

		else:
			prev_marker = marker
			prev_count = prev_count + 1	
		

	if prev_marker == "A":
		print("Adding rainbow colors")
		colors.extend(cm.rainbow(np.linspace(0.0, 1.0, prev_count)))
	elif prev_marker == "B":
		print("Adding bone colors")
		colors.extend(cm.bone(np.linspace(0.0, 1.0, prev_count)))
	else:
		print("Adding copper colors")
		colors.extend(cm.copper(np.linspace(0.0, 1.0, prev_count)))


	fp.close()

	for i in range(total_count):
		colors[i][3] = 1.0

	return colors

def main():
	
	if len(sys.argv) < 4:
		print("Usage: plot_animation.py <list_of_state_files.txt> <path/to/outfile/outfile_base_name> <total frames in gif> <xyz file with domains>\n", file=sys.stderr)
		sys.exit(1)

	xyzData, minmax = read_in_data(sys.argv[1])
	outfilename = sys.argv[2]
	frameCount = int(sys.argv[3])
	domainFile = sys.argv[4]
	nfr = len(xyzData)
	fps = 10 # Frame per sec
	#print("frames is " + str(nfr))

	numPoints = len(xyzData[0]['x'])
	print("numPoints [" + str(numPoints) + "]")

	#cmap = cm.rainbow(np.linspace(0.0, 1.0, numPoints))

	colors = getColors(domainFile)

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	#sct, = ax.plot([],[],[], "bo-", markersize=2)
	#plot(x, y, color='green', linestyle='dashed', marker='o',
	#markerfacecolor='blue', markersize=12).
	#idxs = range(numPoints)
	#colors = cmap[idxs]

	#rcolors = plt.cm.rainbow(range(numPoints))

	#sct, = ax.plot([], [], [], color='red', linestyle='solid')
	#sct = ax.scatter([], [], [], color=cmap[range(numPoints)])
	sct = ax.scatter([], [], [], color=colors, marker='.')
	#sct, = ax.plot(xyzData[0]['x'], xyzData[0]['y'], xyzData[0]['z'], color=cmap[range(numPoints)], linestyle='solid', marker='.', markersize=2)


	multiplier = round(nfr / frameCount, 0)
	title = ax.set_title('Polymer Frame 0')

	multiplier = int(multiplier)

	if multiplier < 1:
		multiplier = 1


	def updateAni(ifrm):
		print('ifrm is [', ifrm, ']')
		
		sct._offsets3d = (xyzData[ifrm*multiplier]['x'], xyzData[ifrm*multiplier]['y'], xyzData[ifrm*multiplier]['z'])
		title.set_text('Polymer Frame ' + str(ifrm) + ' Timestep ' + str(ifrm*multiplier))		
		#sct.set_3d_properties(xyzData[ifrm]['z'])
		return sct
	
	#ax.set_xlim(minmax['min_x'], minmax['max_x'])
	#ax.set_ylim(minmax['min_y'], minmax['max_y'])
	#ax.set_zlim(minmax['min_z'], minmax['max_z'])
	ax.set_xlim(-30,30)
	ax.set_ylim(-30,30)
	ax.set_zlim(-30,30)
	ani = animation.FuncAnimation(fig,updateAni, frameCount, interval=15)
	#ani = animation.FuncAnimation(fig,updateAni,50, interval=50)
	
	ani.save(outfilename + '.gif', writer='imagemagick', fps = fps)

	
	#Writer = animation.writers['ffmpeg']
	#writer = Writer(fps=fps, metadata=dict(artist='Me'), bitrate=1800)
	#ani.save(outfilename + '.mp4', writer='ffmpeg')	
main()
