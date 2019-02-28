import sys
import os
import subprocess


class UsageError(Exception):
        def __init__(self, err):
                self.err = err
                self.usageStatement = "usage: python3 gradientSimDriver.py <COMMAND> <command-specific options>"
        def usage(self):
                print(self.usageStatement, file=sys.stderr)

        def print_err(self):
                print(self.err + "Exiting!", file=sys.stderr)
                self.usage()
                sys.exit(1)



class TemplateFile:
	
	def __init__(self, filename):
		self.filename = filename
		self.data = dict()

	def readInKeyValue(self):
		fp = open(self.filename)

		for line in fp:
			lineAr = line.rstrip().split('\t')
			key = lineAr[0]
			value = lineAr[1:].join("")
			self.data[key] = value 

		fp.close()

	def check_int(arg):
		return int(arg)

	def check_float(arg):
		return float(arg)
	
	def check_dir(arg):
		if not os.path.isdir(arg):
			subprocess.run(['mkdir', '-p', arg])

		return arg		
	
	def check_file(arg):
		if not os.path.isfile(arg):
			raise ValueError
		
		return arg		

class ModuleUnit:
	def __init__(self):
		pass

	def check_int(arg):
		return int(arg)

	def check_float(arg):
		return float(arg)
	
	def check_marker(arg):
		if (len(arg) != 1):
			raise ValueError
		
		return arg			

	

class Domain(ModuleUnit):
	def __init__(self, data):
		ModuleUnit.__init__()

		try:
			self.marker = data['DOMAIN_MARKER']
			self.markerPerc = data['MARKER_PERC']
			self.start = data['DOMAIN_START']
			self.end = data['DOMAIN_END']
		except KeyError:
			raise
		
		check_params()

	def check_params(self):
		try:
			self.marker = ModuleUnit.check_marker(self.marker)
			self.markerPerc = ModuleUnit.check_float(self.markerPerc)
			self.start = ModuleUnit.check_int(self.start)
			self.end = ModuleUnit.check_int(self.end)
	
		except ValueError:
			raise

class Particle(ModuleUnit):
	def __init__(self,data):
		ModuleUnit.__init__()
		
		try:
			self.marker = data['PARTICLE_MARKER']
			self.mass = data['PARTICLE_MASS']
			
			generateInteractions(data['INTERACTIONS'])

	def generateInteractions(self, data):
		self.interactions = list()
		for marker in data:
			self.interactions.append(ParticleInteraction(marker,data[marker]))

class ParticleInteraction(ModuleUnit):
	def __init__(self, marker2, data):
		ModuleUnit.__init__()
		
		try:
			self.marker2 = marker2
			self.bondType = data['BOND_TYPE']
			self.wiggleDist = data['WIGGLE_DIST']
			self.attractForce = data['ATTRACT_FORCE']
			self.attrE = data['ATTR_E']
			self.tailE = data['TAIL_E']
			self.repE = data['REP_E']
			self.repSigma = data['REP_SIGMA']
			self.attrSigma = data['ATTR_SIGMA']
			self.tailSigma = data['TAIL_SIGMA']

		except KeyError:
			raise
		

		self.check_params()
	

	def check_params(self):
		try:
			self.marker2 = ModuleUnit.check_marker(self.marker2)
			self.wiggleDist = ModuleUnit.check_float(self.wiggleDist)
			self.attractForce = ModuleUnit.check_bool(self.attractForce)
			self.attrE = ModuleUnit.check_float(self.attrE)
			self.repE = ModuleUnit.check_float(self.repE)
			self.repSigma = ModuleUnit.check_float(self.repSigma)
			self.attrSigma = ModuleUnit.check_float(self.attrSigma)
			self.tailSigma = ModuleUnit.check_float(self.tailSigma)
	
		except ValueError:
			raise

class BatchTemplate(TemplateFile):
	def __init__(self, filename):
		TemplateFile.__init__(filename)

class PolymerTemplate(TemplateFile):
	def __init__(self, filename):
		TemplateFile.__init__(filename)
		self.readInKeyValue()

		try:
			self.length = self.data['POLYMER_LEN']
			self.defaultMarker = self.data['DEFAULT_MARKER']
			self.savePath = self.data['SAVE_PATH']
			self.saveFilename = self.data['SAVE_FILENAME']
		
			self.domains = self.data['DOMAINS']

		except KeyError:
			raise

		check_params()

	def check_params(self):
		self.length = TemplateFile.check_int(self.length)
		self.defaultMarker = self.check_marker(self.defaultMarker)
		self.savePath = TemplateFile.check_dir(self.savePath)

	def readInKeyValue():
		fp = open(self.filename)

		inDomain = False
		currDomain = dict()
		for line in fp:
			lineAr = line.rstrip().split('\t')
			key = lineAr[0]
			
			if key.startswith('@'):
				if inDomain == True:
					inDomain = False
					self.data['DOMAINS'].append(Domain(currDomain))
					currDomain = dict()
				else:
					inDomain = True

			elif not inDomain:
				value = lineAr[1]
				self.data[key] = value 
			elif inDomain:
				value = lineAr[1]
				currDomain[key] = value 

		fp.close()

	def check_marker(arg):
		if (len(arg) != 1):
			raise ValueError
		
		return arg			
		

class SimTemplate(TemplateFile):
	def __init__(self, filename):
		TemplateFile.__init__(filename)
		self.readInKeyValue()

		try:
			self.energyMinBool = self.data['ENERGY_MIN_BOOL']
			self.energyMinTol = self.data['ENERGY_MIN_TOL']
			self.energyMinMaxIter = self.data['ENERGY_MIN_MAX_ITER']

			self.numBlocks = self.data['NUM_BLOCKS']		
			self.stepsPerBlock = self.data['STEPS_PER_BLOCK']
			
			self.center = self.data['CENTER']
			self.confinementDensity = self.data['CONFINEMENT_DENSITY']
			self.confinementK = self.data['CONFINEMENT_K']
			
			self.stiffness = self.data['STIFFNESS'] 

			self.integrator = self.data['INTEGRATOR']
			self.errorTol = self.data['ERROR_TOL']
			self.thermostat = self.data['THERMOSTAT']
			self.temperature = self.data['TEMPERATURE']

			self.inputPolymerFile = self.data['INPUT_POLYMER_FILE']
			self.savePath = self.data['SAVE_PATH']
			self.saveFilename = self.data['SAVE_FILENAME']
			self.platform = self.data['PLATFORM']

			self.particleTypes = self.data['PARTICLE_TYPE_LIST']

			self.particles = self.data['PARTICLES']

		except KeyError:
			raise

		check_params()

	def check_params(self):

		try:
			self.energyMinBool = TemplateFile.check_bool(self.energyMinBool)
			self.energyMinTol = TemplateFile.check_float(self.energyMinTol)
			self.energyMinMaxIter = TemplateFile.check_int(self.energyMinMaxIter)
		
			self.numBlocks = TemplateFile.check_int(self.numBlocks)
			self.stepsPerBlock = TemplateFile.check_int(self.numBlocks)
		
			self.center = TemplateFile.check_bool(self.center)
			self.confinementDensity = TemplateFile.check_float(self.confinementDensity)
			self.confinementK = TemplateFile.check_int(self.confinementK)
	
			self.stiffness = TemplateFile.check_int(self.stiffness)
		
			self.errorTol = TemplateFile.check_float(self.errorTol)
			self.thermostat = TemplateFile.check_float(self.thermostat)
			self.temperature = TemplateFile.check_float(self.temperature)
		
			self.savePath = TemplateFile.check_dir(self.savePath)
		
			for p in self.particleTypes:
				p = self.check_marker(p)
		except ValueError:
			raise

	def check_marker(arg):
		if (len(arg) != 1):
			raise ValueError
		
		return arg			

	def readInKeyValue():
		fp = open(self.filename)

		inParticle = False
		currParticle = dict()
		for line in fp:
			lineAr = line.rstrip().split('\t')
			key = lineAr[0]
			
			if key.startswith('@'):
				if inParticle == True:
					inParticle = False
					self.data['PARTICLES'].append(Particle(currParticle))
					currParticle = dict()
				else:
					inParticle = True
					marker = lineAr[1]
					currParticle['PARTICLE_MARKER'] = marker

			elif not inParticle:
				value = lineAr[1]
				if key == 'PARTICLE_TYPE_LIST':
					particleTypes = value.split(',')
					self.data['PARTICLE_TYPE_LIST'] = particleTypes

				else:
					self.data[key] = value 
			elif inParticle:
				if key == 'PARTICLE_MASS':
					value = lineAr[1]
					currParticle[key] = value
				else:
					particle2 = lineAr[1]
					value = lineAr[2]
					currParticle['INTERACTIONS'][particle2] = dict()
					currParticle['INTERACTIONS'][particle2][key] = value 

		fp.close()

class RunTemplate(TemplateFile):
	def __init__(self, filename):
		TemplateFile.__init__(filename)
		TemplateFile.readInKeyValue()

		try:
			self.reps = self.data['REPS']
			self.jobName = self.data['SAVE_JOBNAME']
			self.savePath = self.data['SAVE_PATH']
			self.polymerTemplateFile = self.data['POLYMER_TEMPLATE']
			self.simTemplateFile = self.data['SIM_TEMPLATE']
			self.batchTemplateFile = self.data['BATCH_TEMPLATE']
		except KeyError:
			raise

		check_params()

	def check_params(self):
		try:
			self.reps = TemplateFile.check_int(self.reps)
			self.savePath = TemplateFile.check_dir(self.savePath)
			self.polymerTemplateFile = TemplateFile.check_file(self.polymerTemplateFile)
			self.simTemplateFile = TemplateFile.check_file(self.simTemplateFile)
			self.batchTemplateFile = TemplateFile.check_file(self.batchTemplateFile)
		except ValueError:
			raise

	def load_templates(self):
		self.batchParams = BatchTemplate(self.batchTemplateFile)
		self.polymerParams = PolymerTemplate(self.batchTemplateFile)
		self.simParams = SimTemplate(self.batchTemplateFile)

class RunTemplateGradient(RunTemplate):
	def __init__(self, filename):
		RunTemplate.__init__(filename)

		try:
			self.gradientParam = self.data['GRADIENT_PARAM']
			self.gradientStart = self.data['GRADIENT_START']
			self.gradientEnd = self.data['GRADIENT_END']
			self.gradientPoints = self.data['GRADIENT_POINTS']
		except KeyError:
			raise


	def check_params(self):
		RunTemplate.check_params()
		try:
			self.gradientStart = TemplateFile.check_float(self.gradientStart)
			self.gradientEnd = TemplateFile.check_float(self.gradientEnd)
			self.gradientPoints = TemplateFile.check_int(self.gradientPoints)
		except ValueError:
			raise

class SimulationRunGenerator:
	def factory(runType, *args):
		if runType "GRADIENT" and len(args) == 1:
			return RunTemplateGradient(args[0])
		else:
			raise UsageError("Unimplemented feature: runType [" + runType + "] with [" + str(len(args)) + "] args!" )		

	factory = staticmethod(factory)	

	
def main()
	try:
		paramObj = SimulationRunGenerator.factory(sys.argv[1], sys.argv[2:])
	except ValueError:
		UsageError("Something is wrong with your inputs!\n")

	

main()
