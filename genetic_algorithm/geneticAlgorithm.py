'''
Genetic algorithm
'''

import sys
import time
import datetime
import numpy as np
from constants import const
from chromosome import Chromosome
from generation import Generation
from fitnessFunction import get_fit
from graphics.plot_fit import plot_fit, update_fit_plot, close_fit_plot
from graphics.plot_fitness import plot_fitness, update_fitness_plot, close_fitness_plot

def display_generation(generation):
	for k in range(generation.size):
		sys.stdout.write('Chromosome %3d: ' % k)
		for gene in generation.chromosomes[k].genes:
			sys.stdout.write('%6.2f' % gene)
		sys.stdout.write(' rmsd = %6.4f' % generation.chromosomes[k].fitness)
		sys.stdout.write('\n')
	sys.stdout.write('\n')

def estimate_errors(errors, variables, data, threshold):
	# Set the maximal acceptable RMSD
	max_fitness = threshold * data['fitness'][0]
	# Select the sets of the parameters which yeild the RMSD value below the RMSD threshold
	idx = min(range(len(data['fitness'])), key=lambda i: abs(data['fitness'][i]-max_fitness))
	for name in variables:
		var_min = np.amin(data[name][0:idx])
		var_max = np.amax(data[name][0:idx])
		var_dev = 0.5 * (var_max - var_min)
		errors[name] = var_dev

def order_variables(genes, indices, fixed, errors):
	values = []
	optimized = []
	precision = []
	for name in const['variableNames']:
		index = indices[name]
		if not (index == -1):
			values.append(genes[index] / const['variableScales'][name])
			optimized.append('Y')
			if name in errors:
				precision.append(errors[name] / const['variableScales'][name])
			else:
				precision.append(np.nan)
		else:
			values.append(fixed[name] / const['variableScales'][name])
			optimized.append('N')
			precision.append(np.nan)
	parameters = []
	for i in range(len(indices)):
		parameters.append([const['longVariableNames'][i], values[i], optimized[i], precision[i]])
	return parameters	
	
class GeneticAlgorithm:

	def __init__(self, settings, exp, fitPar):
		self.num_generations = settings['num_generations']
		self.size_generation = settings['size_generation']
		self.prob_crossover = settings['prob_crossover']
		self.prob_mutation = settings['prob_mutation']
		self.display_graphics = settings['display_graphics']
		self.best_fit = np.zeros(exp['spc'].size)
		self.best_fitness = np.zeros(self.num_generations)
		self.best_parameters = []
		self.fitness_vs_genes = []
		self.gene_errors = {}

	def run_optimization(self, fitPar, sim, exp, spinA, spinB, calc):
		# Display status
		sys.stdout.write('Starting fitting...\n')
		time_start = time.time()
		# Increment over generations
		for i in range(self.num_generations):
			if (i == 0):
				# Create the first generation
				generation = Generation(self.size_generation)
				generation.first_generation(fitPar['variables']['bounds'])
			else:
				# Create the next generation
				generation.produce_offspring(fitPar['variables']['bounds'], self.prob_crossover, self.prob_mutation)
			# Score the generation
			generation.score_chromosomes(fitPar['variables']['indices'], fitPar['variables']['fixed'], sim, exp, spinA, spinB, calc)
			# Sort chromosomes according to their fitness
			generation.sort_chromosomes() 
			# Save the best fitness in each optimization step
			self.best_fitness[i] = generation.chromosomes[0].fitness
			if self.display_graphics:
				# Determine the best genes and the best fit so far
				best_genes = generation.chromosomes[0].genes
				best_fit = get_fit(best_genes, fitPar['variables']['indices'], fitPar['variables']['fixed'], sim, exp, spinA, spinB)
				# Plot the experimental and simulated spectra
				if (i == 0):
					fig_fit, graph_fit = plot_fit(exp['f'], best_fit, exp['f'], exp['spc'], calc)
					fig_fitness, axes_fitness = plot_fitness(self.best_fitness)
				else:
					update_fit_plot(fig_fit, graph_fit, best_fit)
					update_fitness_plot(axes_fitness, self.best_fitness)
			# Output the current status
			sys.stdout.write('\r')
			sys.stdout.write("Optimization step %d / %d: RMSD = %f" % (i+1, self.num_generations, self.best_fitness[i]))
			sys.stdout.flush()
		if self.display_graphics:
			close_fit_plot(fig_fit)
			close_fitness_plot(fig_fitness)
		# Determine the best genes and the best fit
		best_genes = generation.chromosomes[0].genes	
		self.best_fit = get_fit(best_genes, fitPar['variables']['indices'], fitPar['variables']['fixed'], sim, exp, spinA, spinB)
		# Display status
		time_finish = time.time()
		time_elapsed = str(datetime.timedelta(seconds = time_finish - time_start))
		sys.stdout.write('\n')
		sys.stdout.write('Optimization is finished. Total duration: %s' % (time_elapsed))
		sys.stdout.write('\n\n')
		# Estimate the precision of the genes obtained
		sys.stdout.write('Estimating the precision of fitting parameters... ')
		self.fitness_vs_genes, gene_errors = self.calculate_fitness_vs_genes(best_genes, fitPar, sim, exp, spinA, spinB, calc)
		print gene_errors
		sys.stdout.write('[DONE]\n\n')
		# Make a summary of all fitting parameters
		self.best_parameters = order_variables(best_genes, fitPar['variables']['indices'], fitPar['variables']['fixed'], gene_errors)
		sys.stdout.write('Optimized parameters:\n')
		sys.stdout.write("{0:<16s} {1:<16s} {2:<16s} {3:<16s}\n".format('Parameter', 'Value', 'Optimized', 'Precision (+/-)'))
		for entry in self.best_parameters:
			sys.stdout.write("{0:<16s} {1:<16.2f} {2:<16s} {3:<16.2f}\n".format(entry[0], entry[1], entry[2], entry[3]))
		sys.stdout.write('\n')
	
	def calculate_fitness_vs_genes(self, best_genes, fitPar, sim, exp, spinA, spinB, calc):
		# Initialize arrays for error estimates		
		fitness_vs_genes = []
		gene_errors = {}
		# Set the number of error estimations
		M = len(fitPar['errors']['variables'])
		# Set the number of Monte-Carlo samples in each of calculation
		N = fitPar['errors']['Ns']
		# Increment over different sets of variables
		for i in range(M):
			# Number of variables 
			dim = len(fitPar['errors']['variables'][i])
			# Create a generation
			generation = Generation(N)
			generation.first_generation(fitPar['variables']['bounds'])
			# Set all genes to the optimized values except the ones given in varError  
			for name in const['variableNames']:
				index = fitPar['variables']['indices'][name]
				if not (index == -1) and not (name in fitPar['errors']['variables'][i]):
					for j in range(N):
						generation.chromosomes[j].genes[index] = best_genes[index]
			# Score the generation
			generation.score_chromosomes(fitPar['variables']['indices'], fitPar['variables']['fixed'], sim, exp, spinA, spinB, calc)
			# Sort chromosomes according to their fitness
			generation.sort_chromosomes()
			# Store the variables and the corresponding fitness values
			data = {}
			for name in fitPar['errors']['variables'][i]:
				data[name] = []
				index = fitPar['variables']['indices'][name]
				for j in range(N):
					data[name].append(generation.chromosomes[j].genes[index])
			data['fitness'] = []
			for j in range(N):
				data['fitness'].append(generation.chromosomes[j].fitness)
			fitness_vs_genes.append(data)
			# Estimate the errors
			estimate_errors(gene_errors, fitPar['errors']['variables'][i], data, fitPar['errors']['threshold']) 			
		return [fitness_vs_genes, gene_errors]