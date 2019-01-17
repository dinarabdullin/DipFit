'''
Genetic Algorithm: Plot the results of the fitting
'''

import sys
import numpy as np
from plot_fit import plot_fit
from plot_fitness import plot_fitness
from plot_fitness_vs_genes import plot_fitness_vs_genes
from matplotlib import rcParams

def plot_fitting_data(fit, exp, fitPar, calc, output):
	sys.stdout.write('Plotting the fitting results... ')
	rcParams['font.size'] = 18
	# Plot the fit
	filename = output['directory'] + 'fit.png'
	fig_fit, graph_fit = plot_fit(exp['f'], fit.best_fit, exp['f'], exp['spc'], calc, output['save_figures'], filename)
	# Plot the fitness vs the number of the generation
	filename = output['directory'] + 'fitness.png'
	fig_fitness, axes_fitness = plot_fitness(fit.best_fitness, output['save_figures'], filename)
	# Plot the fitness vs genes
	filename = output['directory'] + 'parameter_errors.png'
	plot_fitness_vs_genes(fit.fitness_vs_genes, fitPar['errors']['variables'], fitPar['variables']['indices'], fitPar['variables']['bounds'], True, filename)	
	sys.stdout.write('[DONE]\n\n')