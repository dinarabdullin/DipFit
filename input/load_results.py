'''
Reads the input data from the config file
'''

import os
import sys
sys.path.append('../')
import argparse
import matplotlib
import numpy as np
from constants import const
from genetic_algorithm.geneticAlgorithm import estimate_errors
from graphics.plot_fit import plot_fit
from graphics.plot_fitness import plot_fitness
from graphics.plot_fitness_vs_genes import plot_fitness_vs_genes
from graphics.additional import keep_figures_live

datasets = [
    {'filename': 'parameter_errors-r_mean.dat',             'variables': ['r_mean']               },
    {'filename': 'parameter_errors-r_width.dat',            'variables': ['r_width']              },
    {'filename': 'parameter_errors-r_mean-r_width.dat',     'variables': ['r_mean', 'r_width']    },
    {'filename': 'parameter_errors-xi_mean.dat',            'variables': ['xi_mean']              },
    {'filename': 'parameter_errors-xi_width.dat',           'variables': ['xi_width']             },
    {'filename': 'parameter_errors-xi_mean-xi_width.dat',   'variables': ['xi_mean', 'xi_width']  },
    {'filename': 'parameter_errors-phi_mean.dat',           'variables': ['phi_mean']             },
    {'filename': 'parameter_errors-phi_width.dat.dat',      'variables': ['phi_width']            },
    {'filename': 'parameter_errors-phi_mean-phi_width.dat', 'variables': ['phi_mean', 'phi_width']},
    {'filename': 'parameter_errors-temp.dat',               'variables': ['temp']                 }
]

def load_fit(path):
	f = []
	spc = []
	fit = []
	filename = path + 'fit.dat'
	if os.path.isfile(filename):
		data = np.genfromtxt(filename, skip_header=1)
		nx = data.shape[0]
		f = [data[j][0] for j in range(nx)]
		spc = [data[j][1] for j in range(nx)]
		fit = [data[j][2] for j in range(nx)]
	return [f, spc, fit]	
		
def load_best_parameters(path):	
	best_parameters = []
	names = []
	values = []
	optimized = []
	precision = []
	filename = path + 'parameters.dat'
	if os.path.isfile(filename):
		file = open(filename, 'r')
		count = 0
		for line in file:
			if (count > 0):
				names.append(str(line[0:15].strip()))
				values.append(float(line[16:31].strip()))
				optimized.append(str(line[32:47].strip()))
				precision.append(float(line[48:63].strip()))
			count = count + 1	
	for i in range(len(names)):
		best_parameters.append([names[i], values[i], optimized[i], precision[i]])
	return best_parameters

def load_fitness_vs_step(path):
	fitness_vs_step = []
	filename = path + 'fitness.dat'
	if os.path.isfile(filename):
		data = np.genfromtxt(filename, skip_header=1)
		nx = data.shape[0]
		fitness_vs_step = [data[j][1] for j in range(nx)]
	return fitness_vs_step
	
def load_fitness_vs_variables(path):
	nf = len(datasets)
	variables = []
	fitness_vs_variables = []
	for i in range(nf):
		filename = path + datasets[i]['filename']
		# Check if the file exists
		if os.path.isfile(filename):
			# If the filr exists, save its ID
			variables.append(datasets[i]['variables'])
			# Read the data from the file
			data = np.genfromtxt(filename, skip_header=1)
			nx = data.shape[0]
			data_dict = {}
			if len(datasets[i]['variables']) == 1:
				name = datasets[i]['variables'][0]
				data_dict[name] = [data[j][0] for j in range(nx)]
				data_dict['fitness'] = [data[j][1] for j in range(nx)]
			elif len(datasets[i]['variables']) == 2:
				name1 = datasets[i]['variables'][0]
				name2 = datasets[i]['variables'][1]
				data_dict[name1] = [data[j][0] for j in range(nx)]
				data_dict[name2] = [data[j][1] for j in range(nx)]
				data_dict['fitness'] = [data[j][2] for j in range(nx)]
			fitness_vs_variables.append(data_dict)
	return [variables, fitness_vs_variables]
	
def compute_errors(variables, fitness_vs_variables):
	parameter_errors = {}
	threshold = 1.1
	M = len(variables)
	for i in range(M):
		estimate_errors(parameter_errors, variables[i], fitness_vs_variables[i], threshold)
	return parameter_errors

if __name__ == '__main__':
	########################
	# Load the DipFit data #
	########################
	# Input arguments
	parser = argparse.ArgumentParser()
	parser.add_argument('path', help="the path to the folder with data files")
	parser.add_argument('-fs', '--fontsize', type=int, help="Set the font size of the figure")
	args = parser.parse_args()	
	# Set the data file name
	path = args.path
	# Set the fontsize
	if args.fontsize:
		fontsize = args.fontsize
		matplotlib.rcParams.update({'font.size': fontsize})
	else:
		matplotlib.rcParams.update({'font.size': 18})
	
	# Load the fit
	[f, spc, best_fit] = load_fit(path)
	print len(f)
	print len(spc)
	print len(best_fit)
	
	# Load the best fitting parameters
	best_parameters = load_best_parameters(path)
	
	# Load fitness vs optimization step
	fitness_vs_step = load_fitness_vs_step(path)
	
	# Load fitness vs individual fitting parameters
	[variables, fitness_vs_variables] = load_fitness_vs_variables(path)
	
	## Recalculate the parameter errors
	# parameter_errors = compute_errors(variables, fitness_vs_variables)
	# print parameter_errors
	
	############################################
	# Make a summary of all fitting parameters #
	############################################
	sys.stdout.write('\nOptimized parameters:\n')
	sys.stdout.write("{0:<16s} {1:<16s} {2:<16s} {3:<16s}\n".format('Parameter', 'Value', 'Optimized', 'Precision (+/-)'))
	for entry in best_parameters:
		sys.stdout.write("{0:<16s} {1:<16.2f} {2:<16s} {3:<16.2f}\n".format(entry[0], entry[1], entry[2], entry[3]))
	sys.stdout.write('\n')
	
	########################
	# Plot the DipFit data #
	########################
	sys.stdout.write('Plotting the fitting results... ')
	# Plot the fit
	plot_fit(f, best_fit, f, spc)
	# Plot the fitness vs the number of the generation
	plot_fitness(fitness_vs_step)
	# Plot the fitness vs genes
	plot_fitness_vs_genes(fitness_vs_variables, variables)	
	sys.stdout.write('[DONE]\n\n')
	# Keep all figures live
	keep_figures_live()
	