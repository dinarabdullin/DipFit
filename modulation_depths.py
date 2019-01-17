'''
Calculates the modulation depth vs g-factor and temperature
'''

import sys
import argparse
import numpy as np
from input import read_config
from spin_func.flipProbabilities import flipProbabilities
from output.output_directory import output_directory
from output.save_flip_probabilities import save_flip_probabilities
from graphics.plot_flip_probabilities import plot_flip_probabilities
from graphics.additional import keep_figures_live
from simulation.simulation import Simulator

def modulation_depth(spin, exp, T):
	# Set the g-axis
	gmin = np.amin(spin['g'])
	gmax = np.amax(spin['g'])
	ginc = 0.001
	Ng = int((gmax - gmin) / ginc) + 1
	g = np.linspace(gmin, gmax, Ng)
	# Calculate the modulation depth vs g and T
	l = []
	for i in range(len(T)):
		p = flipProbabilities(spin['g'], g, exp['magnField'], T[i])
		l.append(p)
	return [g, l]

def mean_modulation_depth(spinA, spinB, exp, calc, T):
	# Spectral simulator
	sim = Simulator(calc)
	# Calculate the effective g-factors
	gA, gB, qA, qB = sim.precalculations(spinA, spinB)
	# Calculate the mean modulation depth vs T
	mean_lambdas = []
	for i in range(len(T)):
		p = flipProbabilities(spinB['g'], gB, exp['magnField'], T[i])
		mean_lambda = np.mean(p)
		mean_lambdas.append(mean_lambda)
	return mean_lambdas

if __name__ == '__main__':
	# Read out the path to the config file
	parser = argparse.ArgumentParser()
	parser.add_argument('filepath', help="The path to the configuration file")
	args = parser.parse_args()
	configPath = args.filepath

	# Read out the config file
	exp, spinA, spinB, mode, simPar, fitPar, calc, output = read_config(configPath)
	# Set the temperature values in K
	T = [3.0, 300.0]

	# Make an output directory
	output_directory(output, configPath)

	# Calculate the modulation depth vs g and T
	sys.stdout.write('Calculating the modulation depth vs g and T... ')
	[g, l] = modulation_depth(spinB, exp, T)
	sys.stdout.write('[DONE]\n\n')

	# Save the data
	if output['save_data']:
		sys.stdout.write('Saving the simulation results into the directory:\n')
		sys.stdout.write(output['directory'])
		filename = output['directory'] + 'modulation_depth.dat'
		save_flip_probabilities(g, l, T, filename)
		sys.stdout.write(' [DONE]\n\n')

	# Plot the data
	sys.stdout.write('Plotting the simulation results... ')
	filename = output['directory'] + 'modulation_depth.png'
	plot_flip_probabilities(g, l, T, output['save_figures'], filename)
	sys.stdout.write('[DONE]\n\n')

	# Calculate the mean modulation depth vs T 
	sys.stdout.write('Calculating the mean modulation depth vs T... ')
	mean_lambdas = mean_modulation_depth(spinA, spinB, exp, calc, T)
	sys.stdout.write('[DONE]\n')
	for i in range(len(T)):
		sys.stdout.write("Mean modulation depth at %s: %f\n" % (str(T[i]) + ' K', mean_lambdas[i]))
	sys.stdout.write('\n')

	# Keep all figures live
	keep_figures_live()