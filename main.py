'''
Simulation of dipolar spectra and dipolar time traces
Requirements:
- python 2.7.X (including modules time, sys, argparse, math, os, datetime, errno, shutil)
- numpy
- scipy
- matplotlib
- libconf
'''

import argparse
from input.input import read_config
from simulation.simulation import Simulator
from genetic_algorithm.geneticAlgorithm import GeneticAlgorithm
from output.output_directory import output_directory
from output.save_simulation_data import save_simulation_data
from output.save_fitting_data import save_fitting_data
from graphics.plot_simulation_data import plot_simulation_data
from graphics.plot_fitting_data import plot_fitting_data
from graphics.additional import keep_figures_live

if __name__ == '__main__':
	# Read out the path to the config file
	parser = argparse.ArgumentParser()
	parser.add_argument('filepath', help="The path to the configuration file")
	args = parser.parse_args()
	configPath = args.filepath

	# Read out the config file 
	exp, spinA, spinB, mode, simPar, fitPar, calc, output = read_config(configPath)

	# Make an output directory
	output_directory(output, configPath)

	# Spectral simulator
	sim = Simulator(calc)

	# Simulation
	if mode['simulation']:
		# Run the simulation
		sim.run_simulation(simPar, exp, spinA, spinB, calc)	
		# Save simulation results
		save_simulation_data(sim, exp, simPar['settings'], output)	
		# Plot simulation results 
		plot_simulation_data(sim, exp, simPar['settings'], calc, output)

	# Fitting
	if mode['fitting']:
		# Init the fitting mode of the simulator
		sim.init_fitting(fitPar, exp, spinA, spinB)
		# Optimizer
		fit = GeneticAlgorithm(fitPar['settings'], exp, fitPar)
		# Run the fitting
		fit.run_optimization(fitPar, sim, exp, spinA, spinB, calc)
		# Plot fitting results
		plot_fitting_data(fit, exp, fitPar, calc, output)
		# Save fitting results
		save_fitting_data(fit, exp, fitPar, output)
		
	# Keep all figures live
	keep_figures_live()