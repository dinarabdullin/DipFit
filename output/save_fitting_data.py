'''
Genetic Algorithm: Saves the results of the fitting
'''

import sys
from save_fit import save_fit
from save_fitness import save_fitness
from save_fitting_parameters import save_fitting_parameters
from save_fitness_vs_genes import save_fitness_vs_genes

def save_fitting_data(fit, exp, fitPar, output):
	if output['save_data']:
		# Output status
		sys.stdout.write('Saving the fitting results into the directory:\n')
		sys.stdout.write(output['directory'])
		# Save the fit to the experimental spectrum
		filename = output['directory'] + 'fit.dat'
		save_fit(exp['f'], fit.best_fit, exp['f'], exp['spc'], filename)
		# Saves the fitness vs optimization step
		filename = output['directory'] + 'fitness.dat'
		save_fitness(fit.best_fitness, filename)
		# Save the optimized fitting parameters
		filename = output['directory'] + 'parameters.dat'
		save_fitting_parameters(fit.best_parameters, filename)
		# Save the fitness vs individual fitting parameters
		save_fitness_vs_genes(fit.fitness_vs_genes, fitPar['errors']['variables'], output['directory'])
		# Output status
		sys.stdout.write(' [DONE]\n\n')  