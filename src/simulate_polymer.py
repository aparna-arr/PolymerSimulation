#!/home/groups/aboettig/Software/anaconda3/envs/openmm-env/bin/python

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import sys
import os
import numpy as np
import subprocess

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
		'INTEGRATOR' : PARAM_OPTS['INTEGRATOR']['variablelangevin'],
		'THERMOSTAT' : 0.02,
		'TEMPERATURE' : 300,
		'INPUT_POLYMER_FILE' : None,
		'SAVE_FILENAME' : "polymer_simulation",
		'SAVE_PATH' : ".",
		'PLATFORM' : PARAM_OPTS['PLATFORM']['cuda'],
		'PARTICLE_TYPE_LIST' : list(), 
		'PARTICLES' : dict(),
		'PARTICLE_MASS' : dict(),
		'REPULSE_FORCE' : PARAM_OPTS['REPULSE_FORCE']['grosberg'],
		'STIFFNESS' : 4,
		'GROSBERG_TRUNC' : 50,
		'ERROR_TOL' : 0.01,
		'ENERGY_MIN_TOL' : 0.3,
		'ENERGY_MIN_MAX_ITER' : 0,
		'NUM_BLOCKS' : 50,
		'STEPS_PER_BLOCK' : 2000
	}
	return params

def init_particle_block(particle, key):
	
	particle = {
		key : {
			'BOND_TYPE' : PARAM_OPTS['BOND_TYPE']['harmonic'],
			'WIGGLE_DIST' : 0.05,
			'ATTRACT_FORCE' : True,
			'ATTR_E' : 0.5,
			'TAIL_E' : 0.1,
			'REP_E' : 3.0,
			'REP_SIGMA' : 1.0,
			'ATTR_SIGMA' : 2.0,
			'TAIL_SIGMA' : 3.0
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
						elif paramKey in ["WIGGLE_DIST", "STIFFNESS", "ATTR_E", "TAIL_E", "REP_E", "REP_SIGMA", "ATTR_SIGMA", "TAIL_SIGMA"]:
							params["PARTICLES"][blockKey][paramKey] = check_float(blockVal)
						elif paramKey in ["ATTRACT_FORCE"]:
							params["PARTICLES"][blockKey][paramKey] = check_bool(blockVal)
						else:
							raise InputError("Input [" + paramKey + "] not an option! Input file [" + filename + "]\n" )		
		
				elif paramKey == "@PARTICLE":
					blockFlag = True
					currParticle = paramVal
				else:
					raise InputError("Input [" + paramKey + "] not an option! Input file [" + filename + "]\n" )		
					
			elif paramKey in ["CONFINEMENT_DENSITY", "CONFINEMENT_K", "THERMOSTAT", "TEMPERATURE", "STIFFNESS", "ERROR_TOL", "ENERGY_MIN_TOL", "GROSBERG_TRUNC"]:
				params[paramKey] = check_float(paramVal)
			elif paramKey in ["ENERGY_MIN_MAX_ITER", "NUM_BLOCKS", "STEPS_PER_BLOCK"]:
				params[paramKey] = check_int(paramVal)
			elif paramKey in ["INTEGRATOR", "PLATFORM", "REPULSE_FORCE"]:
				params[paramKey] = check_param_opt(paramKey,paramVal)
			elif paramKey in ["INPUT_POLYMER_FILE", "SAVE_FILENAME", "SAVE_PATH"]:
				params[paramKey] = paramVal
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

def set_bond_pairs(forceConst):
	# NOTE assumes only ONE chain with linear particle-to-particle bonding!
	for i in range(forceConst['polymerN'] - 1):
		forceConst['bondPairs'].add((i,i+1))

def generate_constant_dictionary(polymerLen, params):
	nm = meter * 1e-9
	conlenScale = 1. # NOTE not sure what conlen is so leaving as constant for now
	kB = BOLTZMANN_CONSTANT_kB * AVOGADRO_CONSTANT_NA

	forceConstants = {
		'nm' : nm,
		'polymerN' : polymerLen,
		'kT' : kB * params['TEMPERATURE'] * kelvin, 
		'conlen' : conlenScale * nm,
		'bondPairs' : set()  
	}	

	set_bond_pairs(forceConstants)

	return forceConstants

def generate_confinement_force(forceConst, params):
	# spherical confinement
	sphereForce = openmm.CustomExternalForce(
		"step(r-SPHaa) * SPHkb * (sqrt((r-SPHaa)*(r-SPHaa) + SPHt*SPHt) - SPHt) "
		";r = sqrt(x^2 + y^2 + z^2 + SPHtt^2)")

	for i in range(forceConst['polymerN']):
		sphereForce.addParticle(i,[])

	radius = (3 * forceConst['polymerN'] / (4 * np.pi * params['CONFINEMENT_DENSITY'])) ** (1 / 3.)

	k = params['CONFINEMENT_K']
	
	sphereForce.addGlobalParameter("SPHkb", k * forceConst['kT'] / forceConst['nm'])
	sphereForce.addGlobalParameter("SPHaa", (radius - 1. / k) * forceConst['nm'])
	sphereForce.addGlobalParameter("SPHt", (1. / k) * forceConst['nm'] / 10.)
	sphereForce.addGlobalParameter("SPHtt", 0.01 * forceConst['nm'])

	return sphereForce

def generate_repulsive_force(forceConst, params):
	# grosberg repulsive force 
	# same between all particles
	repulseOpt = params['REPULSE_FORCE']
	grosbergTrunc = params['GROSBERG_TRUNC']
	radiusGrosberg = forceConst['conlen']
	nbCutOffDist = radiusGrosberg * 2. ** (1. / 6.)

	repul_energy = (
		"step(REPcut2 - REPU) * REPU"
		" + step(REPU - REPcut2) * REPcut2 * (1 + tanh(REPU/REPcut2 - 1));"
		"REPU = 4 * REPe * ((REPsigma/r2)^12 - (REPsigma/r2)^6) + REPe;"
		"r2 = (r^10. + (REPsigma03)^10.)^0.1")
	
	repulseForce = openmm.CustomNonbondedForce(repul_energy)

	repulseForce.addGlobalParameter('REPe', forceConst['kT'])
	repulseForce.addGlobalParameter('REPsigma', radiusGrosberg)
	repulseForce.addGlobalParameter('REPsigma03', 0.3 * radiusGrosberg)

	repulseForce.addGlobalParameter('REPcut', forceConst['kT'] * grosbergTrunc)
	repulseForce.addGlobalParameter('REPcut2', 0.5 * grosbergTrunc * forceConst['kT'])
	
	for i in range(forceConst['polymerN']):
		repulseForce.addParticle(())

	# exclude bonded pairs
	for pair in forceConst['bondPairs']:
		repulseForce.addExclusion(pair[0], pair[1])
		# FIXME should this part be outside the loop?
		repulseForce.setNonbondedMethod(repulseForce.CutoffNonPeriodic)

	return repulseForce

def generate_bond_force(forceConst, params, polymer):
	kbondScalingFactor = float((2 * forceConst['kT'] / (forceConst['conlen']) ** 2) / (kilojoule_per_mole / forceConst['nm'] ** 2))

	# do check here for input... FIXME
	bondForces = openmm.HarmonicBondForce()
	
	for pair in forceConst['bondPairs']:
		pType1 = polymer['type'][pair[0]]
		pType2 = polymer['type'][pair[1]]
	 
		pKey = frozenset([pType1,pType2])

		bondType = params['PARTICLES'][pKey]['BOND_TYPE']
		wiggleDist = params['PARTICLES'][pKey]['WIGGLE_DIST']
		
		# harmonic bonds
		kbond = kbondScalingFactor / (wiggleDist**2)
		bondForces.addBond(pair[0], pair[1], wiggleDist, kbond)	

	return bondForces

def generate_stiffness_force(forceConst, params):
	# add stiffness, same for all particles in polymer
	kStiff = numpy.zeros(forceConst['polymerN'], float) + params['STIFFNESS']
	stiffForce = openmm.CustomAngleForce(
		"kT*angK * (theta - 3.141592) * (theta - 3.141592) * (0.5)")

	for j in range(1, forceConst['polymerN'] - 1):
		stiffForce.addAngle(j - 1, j, j + 1, [float(kStiff[j])])
	
	stiffForce.addGlobalParameter("kT", forceConst['kT'])
	stiffForce.addPerAngleParameter("angK")

	return stiffForce

def construct_attractive_functions(forceConst, polymer, params):
	# generate particle identification string
	particleIdentifyStr = ''
	particleIdentifyFunc = dict()
	particleReverseMap = dict()
	for i in range(len(params['PARTICLE_TYPE_LIST'])):
		particleStr = params['PARTICLE_TYPE_LIST'][i]
		particleIdx = i 
		particleReverseMap[particleStr] = i 
		particleIdentifyStr = particleIdentifyStr + '' + particleStr + ' = ' + str(particleIdx) + ';'
	
		print("particleIdentifyStr [" + str(i) + "] : [" + particleIdentifyStr + "]" )
	
		particleIdentifyFunc[particleStr] = list()
		particleIdentifyFunc[particleStr].append('(step(particleIdx1 - (' + particleStr + '- 0.1)) * step((' + particleStr + '+ 0.1) - particleIdx1))')
		particleIdentifyFunc[particleStr].append('(step(particleIdx2 - (' + particleStr + '- 0.1)) * step((' + particleStr + '+ 0.1) - particleIdx2))')
	
		print("particleStr: [" + particleStr + "]")
		print("particleIdentifyFunc [" + str(i) + "][0] : [" + particleIdentifyFunc[particleStr][0] + "]" )
		print("particleIdentifyFunc [" + str(i) + "][1] : [" + particleIdentifyFunc[particleStr][1] + "]" )
		
	# generate pair identification
	particlePairIdentifyStr = ''
	particlePairDict = dict()
	for i in range(len(params['PARTICLE_TYPE_LIST'])):
		particleType1 = params['PARTICLE_TYPE_LIST'][i];
		particlePairDict[particleType1] = dict()
		for j in range(len(params['PARTICLE_TYPE_LIST'])):
			particleType2 = params['PARTICLE_TYPE_LIST'][j];
			boolPairName = 'is' + particleType1 + particleType2
			particlePairDict[particleType1][particleType2] = boolPairName
			particlePairIdentifyStr = particlePairIdentifyStr + '' + boolPairName + ' = ' + particleIdentifyFunc[particleType1][0] + '*' + particleIdentifyFunc[particleType2][1] + ';'
			
		
	# condense pairs into symmetric pairs (sets)	
	particleCondenseStr = ''
	functionStr = ''
	mainFunctionStr = ''
	condensePairDict = dict()

	# FIXME make this loop readable!
	for i in range(len(params['PARTICLE_TYPE_LIST'])):
		particleType1 = params['PARTICLE_TYPE_LIST'][i];
		for j in range(i, len(params['PARTICLE_TYPE_LIST'])):
			particleType2 = params['PARTICLE_TYPE_LIST'][j];
			boolPair1 = 'is' + particleType1 + particleType2
			boolPair2 = 'is' + particleType2 + particleType1

			pairSymmName = boolPair1 + '_or_' + boolPair2
			particleCondenseStr = particleCondenseStr + '' + pairSymmName + ' = ' + 'step(' + boolPair1 + ' + ' + boolPair2 + ' - 0.1);'

			setKey = frozenset([particleType1,particleType2])
			condensePairDict[setKey] = pairSymmName

			

			ErepulsionPair = 'Erep_' + pairSymmName + ' = rsc12_' + pairSymmName + ' * (rsc2_' + pairSymmName + ' - 1.0) * REPe_' + pairSymmName + ' / emin12 + REPe_' + pairSymmName + ';'

			ErepulsionPair = (ErepulsionPair +\
				'rsc12_' + pairSymmName + ' = rsc4_' + pairSymmName + ' * rsc4_' + pairSymmName + ' * rsc4_' + pairSymmName + ';' +\
				'rsc4_' + pairSymmName + ' = rsc2_' + pairSymmName + ' * rsc2_' + pairSymmName + ';' + \
				'rsc2_' + pairSymmName + ' = rsc_' + pairSymmName + ' * rsc_' + pairSymmName + ';' + \
				'rsc_' + pairSymmName + ' = r / REPsigma_' + pairSymmName + ' * rmin12;')
			
			EattrInnerPair = 'Eattr_inner_' + pairSymmName + ' = - poly_' + pairSymmName + ' * ATTRe_' + pairSymmName + ';'

			EattrOuterPair = 'Eattr_outer_' + pairSymmName + ' = - poly_' + pairSymmName + ' * (ATTRe_' + pairSymmName + ' - TAILe_' + pairSymmName + ') - TAILe_' + pairSymmName + ';'

			polyPair = ('poly_' + pairSymmName + ' = rshft12_' + pairSymmName + ' * (rshft2_' + pairSymmName + ' - 1.0) / emin12 + 1.0' + ';' + \
				'rshft12_' + pairSymmName + ' = rshft4_' + pairSymmName + ' * rshft4_' + pairSymmName + ' * rshft4_' + pairSymmName + ';' +\
				'rshft4_' + pairSymmName + ' = rshft2_' + pairSymmName + ' * rshft2_' + pairSymmName + ';' + \
				'rshft2_' + pairSymmName + ' = rshft_' + pairSymmName + ' * rshft_' + pairSymmName + ';' + \
				'rshft_' + pairSymmName + ' = (r - REPsigma_' + pairSymmName + ' - ATTRdelta_' + pairSymmName + ') / ATTRdelta_' + pairSymmName + ' * rmin12;')

			EtailPair = 'Etail_' + pairSymmName + ' = - TAILe_' + pairSymmName + ' * rtail_' + pairSymmName + ' * rtail_' + pairSymmName + ' * (rtail_' + pairSymmName + ' - 1.0) * (rtail_' + pairSymmName + ' - 1.0) * 16.0;'
		
			EtailPair = EtailPair + 'rtail_' + pairSymmName + ' = (r - REPsigma_' + pairSymmName + ' - 2 * ATTRdelta_' + pairSymmName + ') / TAILr_' + pairSymmName + ' / 2.0 + 0.5;'

			functionStrPair = ErepulsionPair + EattrInnerPair + EattrOuterPair + polyPair + EtailPair
			
			functionStr = functionStr + functionStrPair
 
			mainFuncPair = '(step(REPsigma_' + pairSymmName + ' - r) * Erep_' + pairSymmName +\
				'+ step(r - REPsigma_' + pairSymmName + ') * step(REPsigma_' + pairSymmName + ' + ATTRdelta_' + pairSymmName + ' - r)' \
				' * Eattr_inner_' + pairSymmName + ' * ATTR_bool_' + pairSymmName + \
				'+ step(r - REPsigma_' + pairSymmName + ' - ATTRdelta_' + pairSymmName + ') * step(REPsigma_' + pairSymmName + ' + 2.0 * ATTRdelta_' +\
				pairSymmName + ' - r) * Eattr_outer_' + pairSymmName + ' * ATTR_bool_' + pairSymmName +\
				'+ step(r - REPsigma_' + pairSymmName + ' - ATTRdelta_' + pairSymmName + ') * Etail_' + pairSymmName + ' * ATTR_bool_' + pairSymmName + ') * step(' + pairSymmName + ' - 0.1)'
		
			if mainFunctionStr != "":
				mainFunctionStr = mainFunctionStr + ' + ' + mainFuncPair
			else:
				mainFunctionStr = mainFuncPair
	
	mainFunctionStr += ';'
	equationString = mainFunctionStr + functionStr + particleCondenseStr + particlePairIdentifyStr + particleIdentifyStr 

	print("DEBUG: EQUATION STRING")
	print(equationString)
		
	attractForce = openmm.CustomNonbondedForce(equationString)
	
	# FIXME don't let user set tailRadius
	# as this is used later to set a global parameter (cutoff distance)
	# temporary extremely dumb workaround
	tailRGlobal = 0
	for pair in params['PARTICLES']:
		pairSymmName = condensePairDict[pair]
		attrBool = float(params['PARTICLES'][pair]['ATTRACT_FORCE'])
		attrE = params['PARTICLES'][pair]['ATTR_E']
		tailE = params['PARTICLES'][pair]['TAIL_E']
		repE = params['PARTICLES'][pair]['REP_E']
		repSigma = params['PARTICLES'][pair]['REP_SIGMA']
		attrSigma = params['PARTICLES'][pair]['ATTR_SIGMA']
		tailSigma = params['PARTICLES'][pair]['TAIL_SIGMA']

		print("DEBUG: attrBool: [" + str(attrBool) + "]")
		print("DEBUG: repSigma: [" + str(repSigma) + "]")

		# FIXME continuation of dumb workaround
		tailRGlobal = tailSigma
		
		attractForce.addGlobalParameter('ATTR_bool_' + pairSymmName, attrBool)

		attractForce.addGlobalParameter('REPe_' + pairSymmName, repE * forceConst['kT'])
		attractForce.addGlobalParameter('REPsigma_' + pairSymmName, repSigma * forceConst['conlen'])

		print("Adding param: [" + 'REPsigma_' + pairSymmName + "]")

		attractForce.addGlobalParameter('ATTRe_' + pairSymmName, attrE * forceConst['kT'])
		attractForce.addGlobalParameter('ATTRdelta_' + pairSymmName, forceConst['conlen'] * (attrSigma - repSigma) / 2.0)

		attractForce.addGlobalParameter('TAILr_' + pairSymmName, (tailSigma - attrSigma) * forceConst['kT'])
		attractForce.addGlobalParameter('TAILe_' + pairSymmName, tailE * forceConst['kT'])
		
	# FIXME constants
	attractForce.addGlobalParameter('emin12', 46656.0 / 823543.0)
	attractForce.addGlobalParameter('rmin12', np.sqrt(6.0 / 7.0))
	attractForce.addPerParticleParameter("particleIdx")

	for i in range(len(polymer['type'])):
		currType = polymer['type'][i]
		pidx = particleReverseMap[currType]
		# FIXME something is going wrong with the add Particle thing
		#print("DEBUG[" + str(i) + "]: pidx is [", str(pidx), "] [", str(float(pidx)), "]")
		attractForce.addParticle([pidx])

	# FIXME using the dumb workaround 
	attractForce.setCutoffDistance(forceConst['conlen'] * tailRGlobal)

	# set exclusions for pairs	
	for pair in forceConst['bondPairs']:
		attractForce.addExclusion(pair[0], pair[1])
		attractForce.setNonbondedMethod(attractForce.CutoffNonPeriodic)

	return attractForce

def set_forces_in_system(polymer, params, simSystem):
	# adds all forces to system
	# including polymer forces, bonds, external forces
	# FIXME let user choose grosberg or attractive force

	forceConst = generate_constant_dictionary(len(polymer['type']), params)
	forceList = list()

	forceList.append(generate_confinement_force(forceConst, params))
	# NOTE temporarily commenting out 
	#forceList.append(generate_repulsive_force(forceConst, params))

	forceList.append(construct_attractive_functions(forceConst, polymer, params))

	forceList.append(generate_bond_force(forceConst, params, polymer))
	forceList.append(generate_stiffness_force(forceConst, params))

	for force in forceList:
		simSystem.addForce(force)

def init_system(polymer, params):
	# returns the initialized System
	simSystem = System()
	# add polymer to system
	set_polymer_in_system(polymer, simSystem)
	# add forces to system
	set_forces_in_system(polymer, params, simSystem)

	return simSystem

# integrator functions #
def init_integrator(params):
	ps = second * 1e-12
	# returns the correct integrator object to use with simulation
	if params['INTEGRATOR'] == PARAM_OPTS['INTEGRATOR']['variablelangevin']:
		integ = openmm.VariableLangevinIntegrator(params['TEMPERATURE'] * kelvin, params['THERMOSTAT'] * ( 1/ps ), params['ERROR_TOL'])
		return integ 

# platform functions #
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
	
	# initialize integrator
	integrator = init_integrator(params)

	# initialize platform
	platform = init_platform(params)

	# init simulation
	# using a dummy Topology object because it isn't actually used
	polymerSim = Simulation(Topology(), system, integrator, platform)	

	polymerSim.context.setPositions(polymer['xyz'])

	tolerance = params['ENERGY_MIN_TOL']
	maxIter = params['ENERGY_MIN_MAX_ITER'] 

	# local energy minim
	polymerSim.minimizeEnergy(tolerance, maxIter)

	numBlocks = params['NUM_BLOCKS']
	stepsPerBlock = params['STEPS_PER_BLOCK']

	#checkpoint_str = "_checkpoint_"
	state_str = "_state_"

	if not os.path.isdir(params['SAVE_PATH']):
		subprocess.run(["mkdir", "-p", params['SAVE_PATH']])

	subprocess.run(["cp", params['PARAM_FILE'], params['SAVE_PATH'] + '/'])
	stateListFilename = params['SAVE_PATH'] + '/' + params['SAVE_FILENAME'] + "_state_list.txt"
	statelistFp = open(stateListFilename, 'w')

	for i in range(numBlocks):
		addStr = "block_" + str(i)
		polymerSim.step(numBlocks)
		
		#checkpointFilename = params['SAVE_PATH'] + '/' + params['SAVE_FILENAME'] + checkpoint_str + addStr
		#polymerSim.saveCheckpoint(checkpointFilename)
		stateFilename = params['SAVE_PATH'] + '/' + params['SAVE_FILENAME'] + state_str + addStr
		polymerSim.saveState(stateFilename)

		statelistFp.write(stateFilename + "\n")
	
	statelistFp.close()
main()
