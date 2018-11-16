from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import sys
import os.path
import numpy as np

from enum import Enum
## constants and enums ##
PARAM_OPTS = {
	'INTEGRATOR' : {'variablelangevin' : 1},
	'BOND_TYPE' : {'harmonic' : 1},
	'REPULSE_FORCE' : {'grosberg' : 1},
	'PLATFORM' : {'cuda' : 1}
}

## exception class ##
class UsageError(Exception):

	def __init__(self, err):
		self.err = err
		self.usageStatement = "usage: python3 example_clone.py <simulation_params.txt>"
 
	def usage(self):
		print(self.usageStatement, file=sys.stderr)

	def print_err(self):
		print(self.err + "Exiting!", file=sys.stderr)
		self.usage()
		sys.exit(1)

class InputError(Exception):
	def __init__(self, err):
		self.err = err

	def print_err(self):
		print(self.err + "Exiting!", file=sys.stderr)
		sys.exit(1)

## input functions ##

def check_int(val):
	testInt = int(val)
	return testInt

def check_float(val):
	testFloat = float(val)
	return testFloat

def check_bool(val):
	if (val.lower() == "true" or val.lower() == "t"):	
		return True
	elif (val.lower() == "false" or val.lower() == "f"):
		return False
	else:
		raise ValueError

def make_particle_list(val):
	return [ x.strip() for x in val.split(',') ]

def check_param_opt(opt, val):
	if not opt in PARAM_OPTS:
		raise ValueError
	elif not val.lower() in PARAM_OPTS[opt]:
		raise ValueError

	return PARAM_OPTS[opt][val.lower()]

def check_params(params):
	pass # checks that all required options have been set

def init_params():
	params = {
		'PARAM_FILE' : None,
		'CENTER' : True,
		'CONFINEMENT_DENSITY' : 0.85,
		'CONFINEMENT_K' : 1,
		'ENERGY_MIN' : 10,
		'INTEGRATOR' : PARAM_OPTS['INTEGRATOR']['variablelangevin'],
		'THERMOSTAT' : 0.02,
		'TEMPERATURE' : 300,
		'INPUT_POLYMER_FILE' : None,
		'SAVE_FILENAME' : "polymer_simulation",
		'PLATFORM' : PARAM_OPTS['PLATFORM']['cuda'],
		'PARTICLE_TYPE_LIST' : list(), 
		'PARTICLES' : dict(),
		'PARTICLE_MASS' : dict(),
		'REPULSE_FORCE' : PARAM_OPTS['REPULSE_FORCE']['grosberg'],
		'STIFFNESS' : 4,
		'ERROR_TOL' : 0.01
	}
	return params

def init_particle_block(particle, key):
	
	particle = {
		key : {
			'BOND_TYPE' : PARAM_OPTS['BOND_TYPE']['harmonic'],
			'WIGGLE_DIST' : 0.05
		}
	}
	return particle[key]

def read_in_sim_specs(filename):
	# returns a dictionary of simulation params
	params = init_params()
	params['PARAM_FILE'] = filename

	fp = open(filename, 'r')
	blockFlag = False
	currParticle = "INIT"
	
	try:

		for line in fp:
			if line.startswith('#') or not line.strip():
				continue
	
			#print("LINE IS [" + line + "]")
			paramAr = line.rstrip().split('\t')
			
			# strip comments
			paramKey = paramAr[0].split('#',1)[0].strip()
			paramVal = "\t".join(paramAr[1:]).split('#',1)[0].strip()
			#print("paramKey is [" + paramKey + "]")
			#print("paramVal is [" + paramVal + "]")
	
			if paramKey.startswith('@') or blockFlag == True:
				# beginning of a PARTICLE block
				if paramKey == "@END":
					if paramVal == currParticle:
						blockFlag = False
	
					else:
						raise InputError("End statement [" + paramKey + "] references a different particle than stored as currParticle [" + currParticle + "]. Skipped an @END statement?")				
	
				elif blockFlag == True:
					if paramKey == "PARTICLE_MASS":
						params[paramKey][currParticle] = check_float(paramVal)
					else:
						paramVal = paramVal.split("\t")

						if len(paramVal) != 2:
							raise InputError("Input value for [" + paramKey + "] expected to have 2 column value! Value was [" + " ".join(paramVal) + "] With length [" + str(len(paramVal)) + "] Input file [" + filename + "]\n" )		
	
						otherParticle = paramVal[0].strip()
						blockVal = paramVal[1]
	
						if paramKey in ["BOND_TYPE"]:
							blockKey = frozenset([currParticle, otherParticle])
							params["PARTICLES"][blockKey] = init_particle_block(params["PARTICLES"], blockKey)
							params["PARTICLES"][blockKey][paramKey] = check_param_opt(paramKey, blockVal)
						elif paramKey in ["WIGGLE_DIST", "STIFFNESS"]:
							params["PARTICLES"][blockKey][paramKey] = check_float(blockVal)
						else:
							raise InputError("Input [" + paramKey + "] not an option! Input file [" + filename + "]\n" )		
		
				elif paramKey == "@PARTICLE":
					blockFlag = True
					currParticle = paramVal

					#if params['PARTICLE_TYPE_LIST'] != []:
					#	for i in range(len(params['PARTICLE_TYPE_LIST'])):
					#		if currParticle == params['PARTICLE_TYPE_LIST'][i]:
					#			raise InputError("Error! Two particle blocks for particle [" + currParticle + "]! Input file [" + filename + "]\n")
					#
					#params['PARTICLE_TYPE_LIST'].append(currParticle)
				else:
					raise InputError("Input [" + paramKey + "] not an option! Input file [" + filename + "]\n" )		
					
			elif paramKey in ["CONFINEMENT_DENSITY", "CONFINEMENT_K", "THERMOSTAT", "TEMPERATURE", "STIFFNESS", "ERROR_TOL"]:
				params[paramKey] = check_float(paramVal)
			elif paramKey in ["INTEGRATOR", "PLATFORM", "REPULSE_FORCE"]:
				params[paramKey] = check_param_opt(paramKey,paramVal)
			elif paramKey in ["INPUT_POLYMER_FILE", "SAVE_FILENAME"]:
				params[paramKey] = paramVal
			elif paramKey == "ENERGY_MIN":
				params[paramKey] = check_int(paramVal)
			elif paramKey == "PARTICLE_TYPE_LIST":
				params[paramKey] = make_particle_list(paramVal)
			elif paramKey == "CENTER":
				params[paramKey] = check_bool(paramVal)
			else:
				raise InputError("Input [" + paramKey + "] not an option! Input file [" + filename + "]\n" )		
	
		fp.close()
	
	except Exception as e:
		fp.close()
		raise

	# FIXME: catch all the ValueErrors where they are thrown in this function (in the for loop, in the if/elses
	return params

def read_in_polymer_xyz(filename):
	if not os.path.isfile(filename):
		raise UsageError("Polymer file [" + specsFile + "] does not exist!")
	# returns a Nx3 array of positions, 1 per particle
	# returns a Nx1 array mapping each particle to a particle type
	# final return data structure is a dictionary
	polymer = {
		'xyz' : None,
		'type' : list()
	}

	fp = open(filename, 'r')
	lineCount = 0
	nParticles = 0
	try:
		for line in fp:
			if not line.strip():
				continue
			#print("LINE: [" + line +"]")
			lineAr = line.rstrip().split('\t')
	
			if (len(lineAr) != 4):
				nParticles = int(lineAr[0]) # FIXME error check here
				#print("nParticles is " + str(nParticles))
				polymer['xyz'] = np.zeros((nParticles,3))
				#print("Len is " + str(len(polymer['xyz'])))
				#print("Shape of array is " + str(polymer['xyz'].shape))
				#print("First element is [" + str(polymer['xyz'][0,0]) + "]")
			
				continue
			polymer['type'].append(lineAr[0])

	
			#print("list length " + str(len(polymer['xyz'])))
			# print("lineCount " + str(lineCount))
			polymer['xyz'][lineCount,0] = check_float(lineAr[1]) # x
			polymer['xyz'][lineCount,1] = check_float(lineAr[2]) # y
			polymer['xyz'][lineCount,2] = check_float(lineAr[3]) # z
			lineCount = lineCount + 1
	
		fp.close()
	except Exception:
		fp.close()
		raise

	return polymer

def handle_input(argv):
	if len(argv) != 2:
		raise UsageError("Incorrect number of arguments: [" + str(len(argv) - 1) + "]\n")

	specsFile = argv[1]
	if not os.path.isfile(specsFile):
		raise UsageError("Parameter file [" + specsFile + "] does not exist!")

	try:
		params = read_in_sim_specs(specsFile)
		# polymer is a dictionary
		polymer = read_in_polymer_xyz(params['INPUT_POLYMER_FILE'])	
	except Exception:
		raise
	return polymer, params

## output functions ##
def output_data(filename, data):
	pass

def output_polymer_model(filename, data):
	pass

# polymer functions #
def init_polymer(polymer, params):
	# returns the polymer with centering, randomized (?), and masses updated
	# return is a dictionary

	# set masses
	numParticles = len(polymer['type'])
	polymer['mass'] = list()
	
	for i in range(numParticles):
		particleType = polymer['type'][i]
		particleMass = params['PARTICLE_MASS'][particleType]
		polymer['mass'].append(particleMass)

	# center
	if params['CENTER']:
		avg = np.mean(polymer['xyz'], 0)
		polymer['xyz'] -= avg


	# randomize
	# NOTE skipping this for now. In the openmm-polymer lib, comment says it is to offset coords if integer. I don't know why, so I'm not adding it yet.
	
# system functions #
def set_polymer_in_system(polymer, simSystem):
	# adds polymer positions to system as addParticle
	# add particle masses to system
	for i in range(len(polymer['mass'])):
		simSystem.addParticle(polymer['mass'][i] * 100 * amu)

def set_forces_in_system(polymer, params, simSystem):
	# adds all forces to system
	# including polymer forces, bonds, external forces
	# a.addSphericalConfinement(density=0.85, k=1)
	# a.addHarmonicPolymerBonds(wiggleDist=0.05)
	# a.addGrosbergRepulsiveForce(trunc=50)
	# a.addStiffness(k=4)

	# FIXME constants here!
	polymerN = len(polymer['type'])
	nm = meter * 1e-9
	kB = BOLTZMANN_CONSTANT_kB * AVOGADRO_CONSTANT_NA
	kT = kB * params['TEMPERATURE'] * kelvin 
	conlen = 1. * nm
	bondLen = 1
	kbondScalingFactor = float((2 * kT / (conlen) ** 2) / (kilojoule_per_mole / nm ** 2))
	grosbergTrunc = 50 # FIXME constant! add this to user input

	excludeFromNonbonded = set()
	
	# spherical confinement
	sphereForce = openmm.CustomExternalForce(
		"step(r-SPHaa) * SPHkb * (sqrt((r-SPHaa)*(r-SPHaa) + SPHt*SPHt) - SPHt) "
		";r = sqrt(x^2 + y^2 + z^2 + SPHtt^2)")

	for i in range(polymerN):
		sphereForce.addParticle(i,[])

	radius = (3 * polymerN / (4 * np.pi * params['CONFINEMENT_DENSITY'])) ** (1 / 3.)

	k = params['CONFINEMENT_K']
	
	sphereForce.addGlobalParameter("SPHkb", k * kT / nm)
	sphereForce.addGlobalParameter("SPHaa", (radius - 1. / k) * nm)
	sphereForce.addGlobalParameter("SPHt", (1. / k) * nm / 10.)
	sphereForce.addGlobalParameter("SPHtt", 0.01 * nm)

	# grosberg repulsive force 
	# same between all particles
	repulseOpt = params['REPULSE_FORCE']
	radiusGrosberg = conlen
	nbCutOffDist = radiusGrosberg * 2. ** (1. / 6.)
	repul_energy = (
		"step(REPcut2 - REPU) * REPU"
		" + step(REPU - REPcut2) * REPcut2 * (1 + tanh(REPU/REPcut2 - 1));"
		"REPU = 4 * REPe * ((REPsigma/r2)^12 - (REPsigma/r2)^6) + REPe;"
		"r2 = (r^10. + (REPsigma03)^10.)^0.1")
	
	repulseForce = openmm.CustomNonbondedForce(repul_energy)
	repulseForce.addGlobalParameter('REPe', kT)
	repulseForce.addGlobalParameter('REPsigma', radiusGrosberg)
	repulseForce.addGlobalParameter('REPsigma03', 0.3 * radiusGrosberg)

	repulseForce.addGlobalParameter('REPcut', kT * grosbergTrunc)
	repulseForce.addGlobalParameter('REPcut2', 0.5 * grosbergTrunc * kT)
	
	for i in range(polymerN):
		repulseForce.addParticle(())

	# reference code
        #    for j in range(start, end - 1):
        #        self.addBond(j, j + 1, wiggleDist,
        #            distance=bondLength,
        #            bondType="Harmonic", verbose=False)
        #        if exceptBonds:
        #            self.bondsForException.append((j, j + 1))

	# this part is messy
	# all remaining things require going into the particle struct
	# all-by-all comparisons unless I do it smartly

	# do check here for input... FIXME
	bondForces = openmm.HarmonicBondForce()
	
	for j in range(polymerN - 1):
		pType1 = polymer['type'][j]
		pType2 = polymer['type'][j+1]
	 
		pKey = frozenset([pType1,pType2])

		bondType = params['PARTICLES'][pKey]['BOND_TYPE']
		wiggleDist = params['PARTICLES'][pKey]['WIGGLE_DIST']
		
		# harmonic bonds
		kbond = kbondScalingFactor / (wiggleDist**2)
		bondForces.addBond(j, j+1, wiggleDist, kbond)	
		excludeFromNonbonded.add((j,j+1))
	
	# add stiffness, same for all particles in polymer
	kStiff = numpy.zeros(polymerN, float) + params['STIFFNESS']
	stiffForce = openmm.CustomAngleForce(
		"kT*angK * (theta - 3.141592) * (theta - 3.141592) * (0.5)")

	for j in range(1, polymerN - 1):
		stiffForce.addAngle(j - 1, j, j + 1, [float(kStiff[j])])
	
	stiffForce.addGlobalParameter("kT", kT)
	stiffForce.addPerAngleParameter("angK")

	# check exclusions
	for pair in excludeFromNonbonded:
		#repulseForce.addException(pair[0], pair[1], 0, 0, 0, True)
		repulseForce.addExclusion(pair[0], pair[1])
		repulseForce.setNonbondedMethod(repulseForce.CutoffNonPeriodic)

	simSystem.addForce(sphereForce) # custom external force
	simSystem.addForce(repulseForce) # custom nonbonded force
	simSystem.addForce(bondForces) #exceptBonds = True
	simSystem.addForce(stiffForce) # custom angle force

# FIXME consider setting all the forces in their own functions and taking
# a "factory" style approach to setting them

#def set_stiffness_in_system(params, simSystem):
#	# adds stiffness to system
#	pass

def init_system(polymer, params):
	# returns the initialized System
	simSystem = System()
	# add polymer to system
	set_polymer_in_system(polymer, simSystem)
	# add forces to system
	set_forces_in_system(polymer, params, simSystem)
	# add stiffness to system -> just another force, adding to set_forces
	#set_stiffness_in_system(params, simSystem)

	return simSystem

## # context functions #
## def init_context(system, params)
## 	# initialize context with the platform, integrator, and system
##	pass

# integrator functions #
def init_integrator(params):
	ps = second * 1e-12
	# returns the correct integrator object to use with simulation
	if params['INTEGRATOR'] == PARAM_OPTS['INTEGRATOR']['variablelangevin']:
		integ = openmm.VariableLangevinIntegrator(params['TEMPERATURE'] * kelvin, params['THERMOSTAT'] * ( 1/ps ), params['ERROR_TOL'])
		return integ 

# platform functions
def init_platform(params):
	# returns the correct platform object to use with simulation
	if params['PLATFORM'] == PARAM_OPTS['PLATFORM']['cuda']:
		plat = openmm.Platform.getPlatformByName('CUDA')
		return plat

## MAIN ##

def main():
	try:
		polymer, params = handle_input(sys.argv)
	except UsageError as e:
		e.print_err()
	
	# sets masses, adjusts orientation
	init_polymer(polymer, params)

	# initialize system
	system = init_system(polymer, params)
	
	### initialize context
	##context = init_context(system, params)

	# initialize integrator
	integrator = init_integrator(params)

	# initialize platform
	platform = init_platform(params)

	# init simulation
	# using a dummy Topology object because it isn't actually used
	
	polymerSim = Simulation(Topology(), system, integrator, platform)	

	polymerSim.context.setPositions(polymer['xyz'])

	#FIXME make these constants user options
	tolerance = 0.3
	maxIter = 0 # run until convergence 

	# local energy minim
	polymerSim.minimizeEnergy(tolerance, maxIter)

	# FIXME more constants to make into user options
	numBlocks = 50
	stepsPerBlock = 2000
	checkpoint_str = "_checkpoint_"
	state_str = "_state_"

	for i in range(numBlocks):
		addStr = "block_" + str(i)
		polymerSim.step(numBlocks)
		polymerSim.saveCheckpoint(params['SAVE_FILENAME'] + checkpoint_str + addStr)
		polymerSim.saveState(params['SAVE_FILENAME'] + state_str + addStr)
	
main()

# multiple arrays for each attribute of the polymer
# easy to modify later, smaller data structures ot handle now.
# 1. build polymer

# def read in polymer
## set masses
# def adjust orientation
## center polymer
## randomize orientation
 
# add particle positions to system

# 2. set custom external forces

# def init system
## temperature/thermometer

# def add forces to system
## def add spherical constraints to system
## def add bonds to system
## def add repulsive / attractive forces to system
## def add stiffness to system 

# 3. initialize context
# set platform (CUDA)
# set integrator
# set system

# 4. run simulation
# set the context in sim
# apply local energy minimization
# step through

# 5. save files
# data files
# draw polymer
