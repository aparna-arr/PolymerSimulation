import re
import sys
import os
import subprocess
from enum import Enum

## Constants and Enums ##
PARAM_OPTS = {
	'COMMANDS' : {'ALL' : 1, "GENERATE_BATCH" : 2, "RUN" : 3, "PLOT" : 4}

}

class UsageError(Exception):
        def __init__(self, err):
                self.err = err
                self.usageStatement = "usage: python3 gradientSimDriver.py <COMMAND> <command-specific options>
	Commands:
		ALL
			usage: python3 gradientSimDriver.py GENERATE_TEMPLATES 
				<parameter to gradient over> 
				<gradient start value> 
				<gradient end value> 
				<total number of points> 
				<domain template file> 
				<run template file> 
				<batch template file> 
				<number of repetitions> 
				<output directory> 
				<output prefix>

			-- Parameters to gradient over --
			- any numerical parameter in the domain template 
			- any numerical parameter in the run template
			- POLYMER_LENGTH
			
			For DOMAIN block parameters in the domain template file or PARTICLE block parameters in the run template file, the gradient will apply to ANY parameter which:
				1) matches the name given
				2) has the value $INIT in the template input

		GENERATE_BATCH
	
		RUN
	
		PLOT
"

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


def check_int(val):
	testInt = int(val)
	return testInt

def check_float(val):
	testFloat = float(val)
	return testFloat

def check_file(val):
	if not os.path.isfile(val):
		raise ValueError
	
	return val

def check_dir(val):
	if not os.path.isdir(val):
		subprocess.run(["mkdir", "-p", val])

	return val


def init_params_command_ALL():
	params = {
		'COMMAND' = PARAM_OPTS['COMMANDS']['ALL'],
		'PARAM_TYPE' = 'INIT',
		'PARAM_TO_GRADIENT' = 'INIT',
		'GRADIENT_START' = 'INIT',
		'GRADIENT_END' = 'INIT',
		'GRADIENT_POINTS' = 'INIT',
		'DOMAIN_TEMPLATE' = 'INIT',
		'RUN_TEMPLATE' = 'INIT',
		'BATCH_TEMPLATE' = 'INIT',
		'REP_NUM' = 5,
		'OUTDIR' = '.',
		'OUTPREFIX' = 'GradientSimulation_' 
	}

def _handle_input_all(args):
	if len(args) < 10:
		raise UsageError("Too few arguments for command ALL! [" + str(len(args)) + "] args input.\n")	

	params = init_params_command_ALL()

	try:
		params['PARAM_TO_GRADIENT'] = args[0]
		params['GRADIENT_START'] = check_float(args[1])
		params['GRADIENT_END'] = check_float(args[2])
		params['GRADIENT_POINTS'] = check_int(args[3])
		params['DOMAIN_TEMPLATE'] = check_file(args[4])
		params['RUN_TEMPLATE'] = check_file(args[5])
		params['BATCH_TEMPLATE'] = check_file(args[6])
		params['REP_NUM'] = check_int(args[7])
		params['OUTDIR'] = check_dir(args[8])
		params['OUTPREFIX'] = args[9]
		
		return params

	except Exception as e:
		raise 
	

def handle_input(argv):
	if len(argv) < 2:
		raise UsageError("Too few arguments!\n")

	command = argv[1]
	numArgs = len(argv)
	
	if not command in PARAM_OPTS['COMMANDS']:
		raise UsageError("Command [" + command + "] is not one of the available command options!\n")

	if command == "ALL":
		try: 
			return _handle_input_all(argv[2:])
	elif command == "GENERATE_BATCH":
		raise InputError("Command [" + command + "] is not implemented yet! Sorry!")
	elif command == "RUN":
		raise InputError("Command [" + command + "] is not implemented yet! Sorry!")
	elif command == "PLOT":
		raise InputError("Command [" + command + "] is not implemented yet! Sorry!")
	


def main():
	try:
		params = handle_input(sys.argv)
	except UsageError as e:
		e.print_err()
 
	if params['COMMAND'] == PARAM_OPTS['COMMANDS']['ALL']:
		pass

	## this feels object-oriented ...		

main()

