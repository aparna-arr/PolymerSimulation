import copy
import sys
import os
import subprocess

CHECK_PARAM_TYPES = {
	'DOMAIN' : {
		'DOMAIN_MARKER' : 'marker'
		'MARKER_PERC' : 'float',
		'DOMAIN_START' : 'int',
		'DOMAIN_END' : 'int'
	},
	'PARTICLE' : {
		'PARTICLE_MARKER' : 'marker',
		'PARTICLE_MASS' : 'float'
	},
	'PARTICLE_INTERACTION' : {
		'MARKER' : 'marker',
		'WIGGLE_DIST' : 'float',
		'ATTRACT_FORCE' : 'bool',
		'ATTR_E' : 'float',
		'TAIL_E' : 'float',	
		'REP_E' : 'float',
		'REP_SIGMA' : 'float',
		'ATTR_SIGMA' : 'float',
		'TAIL_SIGMA' : 'float'		
	},
	'POLYMER' : {
		'POLYMER_LEN' : 'int',
		'DEFAULT_MARKER' : 'marker',
		'SAVE_PATH' : 'dir',
	},
	'SIMULATION' : {
		'ENERGY_MIN_BOOL' : 'bool',
		'ENERGY_MIN_TOL' : 'float',
		'ENERGY_MIN_MAX_ITER' : 'int',
		'NUM_BLOCKS' : 'int',
		'STEPS_PER_BLOCK' : 'int',
		'CENTER' : 'bool',
		'CONFINEMENT_DENSITY' : 'float',
		'CONFINEMENT_K' : 'int',
		'STIFFNESS' : 'int',
		'ERROR_TOL' : 'float',
		'THERMOSTAT' : 'float',
		'TEMPERATURE' : 'float'
	},
	'RUN' : {
		'REPS' : 'int',
		'SAVE_PATH' : 'dir',
		'POLYMER_TEMPLATE' : 'file',
		'SIM_TEMPLATE' : 'file',
		'BATCH_TEMPLATE' : 'file'

	},
	'RUN_GRADIENT' : {
		'GRADIENT_START' : 'float',
		'GRADIENT_END' : 'float',
		'GRADIENT_POINTS' : 'int'
	}
}

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
		self.data['SAVE_PATH_TEMPLATE'] = '.'

	def set_savepath(self, savePath):
		try:
			self.data['SAVE_PATH'] = self.check_dir(savePath)
		except ValueError:
			raise

	def set_savepath_template(self, savePath):
		try:
			self.data['SAVE_PATH_TEMPLATE'] = self.check_dir(savePath)
		except ValueError:
			raise

	def set_savefilename(self, saveFilename):
		self.data['SAVE_FILENAME'] = saveFilename

	def set_savefilename_template(self, saveFilename):
		self.data['SAVE_FILENAME_TEMPLATE'] = saveFilename

	def readInKeyValue(self):
		fp = open(self.filename)

		for line in fp:
			lineAr = line.rstrip().split('\t')
			key = lineAr[0]
			value = lineAr[1:].join("")
			self.data[key] = value 

		fp.close()

	def check(self, key):
		for param in CHECK_PARAM_TYPES[key].keys():
			checkType = REQUIRED_PARAMS[key][param]
			
			if checkType == 'marker':
				self.data[param] = self.check_marker(self.data[param])
			elif checkType == 'int':
				self.data[param] = self.check_int(self.data[param])
			elif checkType == 'float':
				self.data[param] = self.check_float(self.data[param])
			elif checkType == 'bool':
				self.data[param] = self.check_bool(self.data[param])
			elif checkType == 'dir':
				self.data[param] = self.check_dir(self.data[param])
			elif checkType == 'file':
				self.data[param] = self.check_file(self.data[param])

	def check_int(self, arg):
		return int(arg)

	def check_float(self, arg):
		return float(arg)
	
	def check_dir(self, arg):
		if not os.path.isdir(arg):
			subprocess.run(['mkdir', '-p', arg])

		return arg		
	
	def check_file(self, arg):
		if not os.path.isfile(arg):
			raise ValueError
		
		return arg		

	def check_marker(self, arg):
		if (len(arg) != 1):
			raise ValueError
		
		return arg			

	def check_bool(self, arg):
		if arg.lower() == 'true':
			return True
		elif arg.lower() == 'false':
			return False
		else:
			raise ValueError

	def write_file(self, contentStr):
		filename = self.data['SAVE_PATH_TEMPLATE'] + '/' + self.data['SAVE_FILENAME_TEMPLATE']
		fp = open(filename, "w")

		fp.write(contentStr)
	
		fp.close()

class ModuleUnit:
	def __init__(self):
		self.data = dict()

	def check(key):
		for param in CHECK_PARAM_TYPES[key].keys():
			checkType = REQUIRED_PARAMS[key][param]
			
			if checkType == 'marker':
				self.data[param] = self.check_marker(self.data[param])
			elif checkType == 'int':
				self.data[param] = self.check_int(self.data[param])
			elif checkType == 'float':
				self.data[param] = self.check_float(self.data[param])
			elif checkType == 'bool':
				self.data[param] = self.check_bool(self.data[param])

	def check_int(arg):
		return int(arg)

	def check_float(arg):
		return float(arg)
	
	def check_marker(arg):
		if (len(arg) != 1):
			raise ValueError
		
		return arg			
	
	def check_bool(arg):
		if arg.lower() == 'true':
			return True
		elif arg.lower() == 'false':
			return False
		else:
			raise ValueError

class Domain(ModuleUnit):
	def __init__(self, data):
		ModuleUnit.__init__()

		self.data = deepcopy(data)	

		try:
			self.check('DOMAIN')	
		except ValueError or KeyError:
			raise

	def generate_string(self):
		domainStr = "@DOMAIN\n"
		
		for param in self.data['KEY_ORDER']:
			domainStr += param + "\t" + str(self.data[param]) + "\n"

		domainStr += "@DOMAIN\n"

		return domainStr

	def modify_param(self, location, param, value):
		if location == 'DOMAIN' and self.data[param] == -1:
			self.data[param] = value

class Particle(ModuleUnit):
	def __init__(self,data):
		ModuleUnit.__init__()
	
		try:
			self.data = deepcopy(data)		
			self.generateInteractions()
			self.check('PARTICLE')
		except ValueError or KeyError:
			raise

	def generateInteractions(self):
		interactions = list()

		for marker in self.data['INTERACTIONS']:
			interactions.append(ParticleInteraction(marker,self.data['INTERACTIONS'][marker], self.data['KEY_ORDER'][marker]))

		self.data['INTERACTIONS'] = deepcopy(interactions)

	def generate_string(self):
		particleStr = "@PARTICLE\t" + self.data['PARTICLE_MARKER'] + "\n"
		particleStr += "PARTICLE_MASS\t" + str(self.data['PARTICLE_MASS']) + "\n"
		for i in self.data['INTERACTIONS']:
			particleStr += self.data['INTERACTIONS'].generate_string() . "\n"

		particleStr += "@END\n"

		return particleStr

	def modify_param(self, location, param, value):
		if location == 'PARTICLE' and self.data[param] == -1:
			self.data[param] = value
		elif location == 'PARTICLE_INTERACTION':
			for i in range(len(self.data['INTERACTIONS'])):
				self.data['INTERACTIONS'][i].modify_param(location, param, value)
					
class ParticleInteraction(ModuleUnit):
	def __init__(self, marker2, data, keyOrder):
		ModuleUnit.__init__()
		
		try:
			self.keyOrder = deepcopy(keyOrder)
			self.data = deepcopy(data)
			self.data['MARKER'] = marker2
			self.check('INTERACTION')
		except ValueError or KeyError:
			raise
	
	def generate_string():
		interactionStr = ""
		marker = self.data['MARKER']
		for param in self.keyOrder:
			interactionStr += param + "\t" + marker + "\t" + str(self.data[param]) + "\n"

		return interactionStr

	def modify_param(self, location, param, value):
		if location == 'PARTICLE_INTERACTION' and self.data[param] == -1:
			self.data[param] = value

class BatchTemplate(TemplateFile):
	def __init__(self, filename):
		TemplateFile.__init__(filename)
		self.data['SAVE_FILENAME_TEMPLATE'] = 'batch.sh'	
		self.readIn()

	def readIn(self):
		fp = open(self.filename)
		
		for line in fp:
			self.data['BATCH_STRING'] += line

		fp.close()

	def set_jobname(self,jobName):
		self.data['BATCH_STRING'] = self.data['BATCH_STRING'].replace('$INIT', jobName)

	def generate_file(self, input_path):
		command = '\nsimulate_polymer ' + input_path + '\n'

		batchStr = self.data['BATCH_STRING'] + command

		self.write_file(batchStr)

class PolymerTemplate(TemplateFile):
	def __init__(self, filename):
		TemplateFile.__init__(filename)
		self.data['SAVE_FILENAME_TEMPLATE'] = 'polymer.polymer'	

		self.keyOrder = list()
		self.readInKeyValue()

		try:
			self.check('POLYMER')
		except ValueError or KeyError:
			raise

	def get_polymer_path(self):
		return self.data['SAVE_PATH'] + '/' + self.data['SAVE_FILENAME']

	def readInKeyValue(self):
		fp = open(self.filename)

		inDomain = False
		currDomain = dict()
		currDomain['KEY_ORDER'] = list()
		for line in fp:
			lineAr = line.rstrip().split('\t')
			key = lineAr[0]
			
			if key.startswith('@'):
				if inDomain == True:
					inDomain = False
					self.data['DOMAINS'].append(Domain(currDomain))
					currDomain = dict()
					currDomain['KEY_ORDER'] = list()
				else:
					inDomain = True

			elif not inDomain:
				value = lineAr[1]
				self.data[key] = value 
				self.keyOrder.append(key)
			elif inDomain:
				value = lineAr[1]
				currDomain[key] = value 
				currDomain['KEY_ORDER'].append(key)

		fp.close()

	def generate_file(self):
		polymerStr = ""
		for param in self.keyOrder:
			polymerStr += param + "\t" + self.data[param] + "\n"
			
		polymerStr += "\n"

		for d in self.data['DOMAINS']:
			polymerStr += d.generate_string() + "\n"

		self.write_file(polymerStr)

	def modify_param(self, location, param, value):
		if location == 'POLYMER':
			self.data[param] = value
		elif location == 'DOMAIN':
			self.data['DOMAINS'].modify_param(location, param, value)

class SimTemplate(TemplateFile):
	def __init__(self, filename):
		TemplateFile.__init__(filename)
		self.keyOrder = list()
		self.data['SAVE_FILENAME_TEMPLATE'] = 'simulation.sim'	

		self.readInKeyValue()

		try:
			self.check_params()
		except ValueError or KeyError:
			raise


	def set_input_polymer(self, polymerPath):
		self.data['INPUT_POLYMER_FILE'] = polymerPath
		
	def check_params(self):
		try:
			self.check('SIMULATION')	
			
			if 'INPUT_POLYMER_FILE' in self.data.keys():
				self.check_file(self.data['INPUT_POLYMER_FILE'])
			
			if 'SAVE_PATH' in self.data.keys():
				self.check_dir(self.data['SAVE_PATH'])

			for p in self.data['PARTICLE_TYPE_LIST']:
				p = self.check_marker(p)

		except ValueError:
			raise

	def readInKeyValue():
		fp = open(self.filename)

		inParticle = False
		currParticle = dict()
		currParticle['INTERACTIONS'] = dict() 
		currParticle['PARTICLE_TYPE_LIST'] = list() 
		currParticle['KEY_ORDER'] = dict()
		for line in fp:
			lineAr = line.rstrip().split('\t')
			key = lineAr[0]
			
			if key.startswith('@'):
				if inParticle == True:
					inParticle = False
					self.data['PARTICLES'].append(Particle(currParticle))
					currParticle = dict()
					currParticle['INTERACTIONS'] = dict() 
					currParticle['PARTICLE_TYPE_LIST'] = list() 
					currParticle['KEY_ORDER'] = dict()
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
					keyOrder.append(key)
			elif inParticle:
				if key == 'PARTICLE_MASS':
					value = lineAr[1]
					currParticle[key] = value
				else:
					particle2 = lineAr[1]
					value = lineAr[2]
					currParticle['INTERACTIONS'][particle2] = dict()
					currParticle['INTERACTIONS'][particle2][key] = value
					currParticle['KEY_ORDER'][particle2].append(key) 

		fp.close()

	def generate_file(self):
		simFileStr = ""

		for param in self.keyOrder :
			if param == 'PARTICLE_TYPE_LIST':
				continue
	
			simFileStr += param + "\t" + self.data[param] + "\n"

		simFileStr += 'PARTICLE_TYPE_LIST\t' + self.data['PARTICLE_TYPE_LIST'].join(',') + '\n'
		simFileStr += '\n'
		
		for p in self.data['PARTICLES'] :
			simFileStr += p.generate_string()

		self.write_file(simFileStr)

	def modify_param(self, location, param, value):
		if location == 'SIMULATION':
			self.data[param] = value
		elif location in ['PARTICLE', 'PARTICLE_INTERACTION']:
			for i in range(len(self.data['PARTICLES'])):
				self.data['PARTICLES'][i].modify_param(location, param, value)

class RunTemplate(TemplateFile):
	def __init__(self, filename):
		TemplateFile.__init__(filename)
		self.readInKeyValue()
		self.data['SAVE_FILENAME_TEMPLATE'] = 'run.run'	

		try:
			self.check('RUN')
		except KeyError:
			raise


		self.batches = list()
		self.polymers = list()
		self.simulations = list()

		self.load_templates()
#		self.set_paths()

	def load_templates(self):
		self.batchParams = BatchTemplate(self.batchTemplateFile)
		self.polymerParams = PolymerTemplate(self.batchTemplateFile)
		self.simParams = SimTemplate(self.batchTemplateFile)

	
#	def set_paths(self):
#		self.batchParams.set_jobname(self.jobName)
#		self.simParams.set_savepath(self.savePath)
#		self.simParams.set_savefilename(self.jobName)

		#self.simParams.set_input_polymer(self.polymerParams.get_polymer_path())
	def generate_run_templates(self):
		for i in range(self.data['REPS']):
			jobName = self.data['SAVE_JOBNAME'] + "_rep_" + str(i)
			batches.append(deepcopy(self.batchParams))
			polymers.append(deepcopy(self.polymerParams))
			simulations.append(deepcopy(self.simParams))			
			# need to do this one first!!
			polymers[i].set_savefilename(jobName + '.xyz')
			polymers[i].set_savepath(self.data['SAVE_PATH'] + '/xyz')

			polymers[i].set_savefilename_template(jobName + '.polymertemplate')
			polymers[i].set_savepath_template(self.data['SAVE_PATH'] + '/templates')

			batches[i].set_jobname(jobName)
			batches[i].set_savefilename(jobName + '.sh')
			batches[i].set_savepath(self.data['SAVE_PATH'] + '/batch')

			simulations[i].set_savefilename(jobName + '.sim')
			simulations[i].set_savepath(self.data['SAVE_PATH'] + '/simulation')
			simulations[i].set_input_polymer(self.polymers[i].get_polymer_path())

			simulations[i].set_savefilename_template(jobName + '.simtemplate')
			simulations[i].set_savepath_template(self.data['SAVE_PATH'] + '/templates')
			 	
	def generate_run(self):
		self.generate_run_templates()
		for i in range(len(self.simulations)):
			self.polymers[i].generate_file()
			self.simulations[i].generate_file()
			self.batches[i].generate_file(self.simulations[i].get_file_path())

	## FIXME write PolymerTemplate.generate_polymer()	
	def generate_polymers(self):
		for i in range(len(self.polymers[i])):
			polymers[i].generate_polymer()

	def generate_master_script(self):
		filename = self.data['SAVE_PATH'] + '/MasterScript_' + self.data['SAVE_JOBNAME'] + '.sh'

		callstr = '#!/bin/bash\n'
		for i in range(len(self.batches)):
			callstr += 'sbatch ' + self.batches[i].get_file_path() . "\n"
		
		fp = open(filename,'w')

		fp.write(callstr)

		fp.close()
	
class RunTemplateGradient(RunTemplate):
	def __init__(self, filename):
		RunTemplate.__init__(filename)

		#self.readInKeyValue()
		try:
			self.check('RUN_GRADIENT')
			self.gradientParamFound = self.find_gradient_param()
		except ValueError or KeyError:
			raise

		self.gradientList = self.generate_gradient()

	def generate_gradient(self):
		step = abs(self.gradientEnd - self.gradientStart) / self.gradientPoints
		
		i = self.gradientStart

		gradientList = list()
		while i <= self.gradientEnd
			gradientList.append(i)
			i += step

		return gradientList
	
	def find_gradient_param(self):

		found = False
		for key in CHECK_PARAM_TYPES.keys():
			if self.gradientParam in CHECK_PARAM_TYPES[key].keys():
				found = key
				
		if not found:
			raise KeyError("Could not find gradient param [" + self.gradientParam + "]")

		return found

	def generate_run_templates(self):
		gradientList = generate_gradient()
		foundParam = find_gradient_param()
		gradParam = self.data['GRADIENT_PARAM']

		for g in gradientList:
			currJobName = self.data['SAVE_JOBNAME'] + "_" + self.data['GRADIENT_PARAM'] + "_" + str(g) + "_"

			currBatch = deepcopy(self.batchParams)
			currPolymer = deepcopy(self.polymerParams)
			currSim = deepcopy(self.simParams)

			if foundParam in ['SIMULATION', 'PARTICLE', 'PARTICLE_INTERACTION']:
				currSim.modify_param(foundParam, gradParam, g)
			elif foundParam in ['POLYMER', 'DOMAIN']:
				currPolymer.modify_param(foundParam, gradParam, g)
						 
			for i in range(self.data['REPS']):
				jobName = currJobName + "_rep_" + str(i)
				batches.append(deepcopy(currBatch))
				polymers.append(deepcopy(currPolymer))
				simulations.append(deepcopy(currSim))			
				# need to do this one first!!
				polymers[i].set_savefilename(jobName + '.xyz')
				polymers[i].set_savepath(self.data['SAVE_PATH'] + '/xyz')

				polymers[i].set_savefilename_template(jobName + '.polymertemplate')
				polymers[i].set_savepath_template(self.data['SAVE_PATH'] + '/templates')

				batches[i].set_jobname(jobName)
				batches[i].set_savefilename(jobName + '.sh')
				batches[i].set_savepath(self.data['SAVE_PATH'] + '/batch')

				simulations[i].set_savefilename(jobName + '.sim')
				simulations[i].set_savepath(self.data['SAVE_PATH'] + '/simulation')
				simulations[i].set_input_polymer(self.polymers[i].get_polymer_path())

				simulations[i].set_savefilename_template(jobName + '.simtemplate')
				simulations[i].set_savepath_template(self.data['SAVE_PATH'] + '/templates')

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

	
	paramObj.prepare()

main()
