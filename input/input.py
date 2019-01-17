'''
Reads the input data from the config file
'''

import io
import sys
import libconf
import numpy as np
from constants import const
from spin_func.gFactorHighSpinFe import gFactorHighSpinFe
from math_func.findMaxOverInterval import findMaxOverInterval

def load_spectrum(filepath):
	xList = []
	yList = []
	file = open(filepath, 'r')
	for line in file:
		str = line.split()
		xList.append(float(str[0]))
		yList.append(float(str[1]))
	xExp = np.array(xList)	
	yExp = np.array(yList)	
	return [xExp, yExp]

def load_timetrace(filepath):
	xList = []
	yList = []
	file = open(filepath, 'r')
	for line in file:
		str = line.split()
		xList.append(float(str[0]))
		yList.append(float(str[1]))
	xExp = np.array(xList)	
	yExp = np.array(yList)	
	return [xExp, yExp]    

def fitting_variables(object):
	indices = {}
	bounds = []
	fixed = {}
	count = 0
	for name in const['variableNames']:
		attribute = getattr(object, name)
		vopt = int(attribute.opt)
		if vopt:
			indices[name] = count
			count += 1
			vmin = float(attribute.range[0]) * const['variableScales'][name]
			vmax = float(attribute.range[1]) * const['variableScales'][name]
			bounds.append([vmin, vmax])
		else:
			indices[name] = -1
			fixed[name] = float(attribute.value) * const['variableScales'][name]
	return [indices, bounds, fixed, count]
	
def read_config(filepath):
	# Display status
	sys.stdout.write('\nReading out the configuration file... ')
	# Initialize dictionaries for input parameters 
	exp = {}
	spinA = {}
	spinB = {}
	mode = {}
	simPar = {}
	fitPar = {}
	calc = {}
	output = {}
	# Read out the config file
	with io.open(filepath) as file:
		config = libconf.load(file)
		###########################
		# Experimental parameters #
		###########################
		exp['path_spectrum'] = config.path_spectrum
		exp['path_timetrace'] = config.path_timetrace
		exp['t'] = []
		exp['sig'] = []
		exp['f'] = []
		exp['spc'] = []
		if exp['path_spectrum']:
			exp['f'], exp['spc'] = load_spectrum(exp['path_spectrum'])
		if exp['path_timetrace']:
			exp['t'], exp['sig'] = load_timetrace(exp['path_timetrace'])
		exp['mwFreq'] = float(config.mw_frequency)
		exp['magnField'] = float(config.magnetic_field)
		##########################
		# Spin system parameters #
		##########################
		spinA['label'] = config.spinA.label
		spinA['D'] = const['wn2MHz'] * np.array([float(config.spinA.D[0]), 
												 float(config.spinA.D[1])])
		spinA['ga'] = np.array([float(config.spinA.g[0]), 
								float(config.spinA.g[1]), 
								float(config.spinA.g[2])])
		spinA['g'] = spinA['ga']
		spinA['gFrame'] = const['deg2rad'] * np.array([float(config.spinA.gFrame[0]), 
													   float(config.spinA.gFrame[1]), 
													   float(config.spinA.gFrame[2])])
		spinB['label'] = config.spinB.label
		spinB['D'] = const['wn2MHz'] * np.array([float(config.spinB.D[0]), 
												 float(config.spinB.D[1])])
		spinB['ga'] = np.array([float(config.spinB.g[0]), 
								float(config.spinB.g[1]), 
								float(config.spinB.g[2])])
		spinB['g'] = spinB['ga']
		if (spinB['label'] == 'hs_fe'):
			spinB['g'] = gFactorHighSpinFe(spinB, exp['mwFreq'])
		spinB['gFrame'] = const['deg2rad'] * np.array([float(config.spinB.gFrame[0]), 
													   float(config.spinB.gFrame[1]), 
													   float(config.spinB.gFrame[2])])
		####################
		# Calculation mode #
		####################
		switch = int(config.mode)
		if (switch == 0):
			mode['simulation'] = 1
			mode['fitting'] = 0
		elif (switch == 1):
			mode['simulation'] = 0
			mode['fitting'] = 1
		###################
		# Simulation mode #
		###################
		if mode['simulation']:
			# Simulation settings
			simPar['settings'] = {}
			simPar['settings']['spc'] = int(config.simulation_settings.spc)
			simPar['settings']['spc_vs_theta'] = int(config.simulation_settings.spc_vs_theta)
			simPar['settings']['spc_vs_xi'] = int(config.simulation_settings.spc_vs_xi)
			simPar['settings']['spc_vs_phi'] = int(config.simulation_settings.spc_vs_phi)
			simPar['settings']['spc_vs_E'] = int(config.simulation_settings.spc_vs_E)
			simPar['settings']['timetrace'] = int(config.simulation_settings.timetrace)
			if not (exp['f'] == []):
				simPar['settings']['faxis_normalized'] = 0
			else:
				simPar['settings']['faxis_normalized'] = int(config.simulation_settings.faxis_normalized)
			# Simulation variables
			simPar['variables'] = {}
			simPar['variables']['r_mean'] = float(config.simulation_variables.r_mean)
			simPar['variables']['r_width'] = float(config.simulation_variables.r_width)
			simPar['variables']['xi_mean'] = const['deg2rad'] * float(config.simulation_variables.xi_mean)
			simPar['variables']['xi_width'] = const['deg2rad'] * float(config.simulation_variables.xi_width)
			simPar['variables']['phi_mean'] = const['deg2rad'] * float(config.simulation_variables.phi_mean)
			simPar['variables']['phi_width'] = const['deg2rad'] * float(config.simulation_variables.phi_width)
			simPar['variables']['temp'] = float(config.simulation_variables.temp)
		################
		# Fitting mode #
		################
		if mode['fitting']:
			# Fitting settings
			fitPar['settings'] = {}
			fitPar['settings']['num_generations'] = int(config.genetic_algorithm.num_generations)
			fitPar['settings']['size_generation'] = int(config.genetic_algorithm.size_generation)
			fitPar['settings']['prob_crossover'] = float(config.genetic_algorithm.prob_crossover)
			fitPar['settings']['prob_mutation'] = float(config.genetic_algorithm.prob_mutation)
			fitPar['settings']['display_graphics'] = int(config.genetic_algorithm.display_graphics)
			# Fitting variables
			fitPar['variables'] = {}
			fitPar['variables']['indices'] = {}
			fitPar['variables']['bounds'] = []
			fitPar['variables']['fixed'] = {}
			fitPar['variables']['indices'], fitPar['variables']['bounds'], fitPar['variables']['fixed'], fitPar['variables']['size'] = fitting_variables(config.fitting_variables)
			# Error estimates
			fitPar['errors'] = {}
			fitPar['errors']['variables'] = config.error_estimation.variables
			fitPar['errors']['Ns'] = int(config.error_estimation.Ns)
			fitPar['errors']['threshold'] = float(config.error_estimation.threshold) * 0.01
		########################
		# Calculation settings #
		########################
		calc['r_distr'] = config.calculation_settings.r_distr
		calc['xi_distr'] = config.calculation_settings.xi_distr
		calc['phi_distr'] = config.calculation_settings.phi_distr
		calc['Ns'] = int(config.calculation_settings.Ns) 
		calc['f_min'] = float(config.calculation_settings.fmin)
		calc['f_max'] = float(config.calculation_settings.fmax)
		if (exp['f'] == []):
			calc['spc_max'] = 1.0
		else:
			if not calc['f_max']:
				calc['f_max'] = np.amax(exp['f'])
				calc['spc_max'] = findMaxOverInterval(exp['spc'],exp['f'],calc['f_min'],calc['f_max'])
			else:
				calc['spc_max'] = findMaxOverInterval(exp['spc'],exp['f'],calc['f_min'],calc['f_max'])
		###############
		# Output data #
		###############
		output['directory'] = config.output.directory
		output['save_data'] = int(config.output.save_data)
		output['save_figures'] = int(config.output.save_figures)
		# Display status
		sys.stdout.write('[DONE]\n')
		if (spinB['label'] == 'hs_fe'):
			sys.stdout.write('Effective g-values of the high-spin Fe(III): gxx = {0:.2f}, gyy = {1:.2f}, gzz = {2:.2f}\n'.format(spinB['g'][0], spinB['g'][1], spinB['g'][2]))
		sys.stdout.write('\n')
	return [exp, spinA, spinB, mode, simPar, fitPar, calc, output]