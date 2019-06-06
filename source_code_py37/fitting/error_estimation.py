'''
Estimate the error of optimized fitting parameters
'''

import numpy as np
from fitting.generation import Generation
from supplement.constants import const

	
def calculate_score_vs_parameters(best_parameters, fitPar, sim, exp, spinA, spinB, calc):
	# Initialize arrays for score vs parameters		
	score_vs_parameters = []
	parameter_errors = {}
	# Set the number of parameter sets
	M = len(fitPar['errors']['variables'])
	# Set the number of Monte-Carlo samples in each calculation
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
					generation.chromosomes[j].genes[index] = best_parameters[name]['value']
		# Score the generation
		generation.score_chromosomes(fitPar['variables']['indices'], fitPar['variables']['fixed'], sim, exp, spinA, spinB, calc)
		# Sort chromosomes according to their score
		generation.sort_chromosomes()
		# Store the variables and the corresponding score values
		score_vs_variables = {}
		for name in fitPar['errors']['variables'][i]:
			score_vs_variables[name] = []
			index = fitPar['variables']['indices'][name]
			for j in range(N):
				score_vs_variables[name].append(generation.chromosomes[j].genes[index])
		score_vs_variables['score'] = []
		for j in range(N):
			score_vs_variables['score'].append(generation.chromosomes[j].score)
		score_vs_parameters.append(score_vs_variables)
		# Estimate the parameter errors
		estimate_parameter_errors(parameter_errors, fitPar['errors']['variables'][i], score_vs_variables, fitPar['errors']['threshold'], best_parameters) 			
	return [score_vs_parameters, parameter_errors]	

	
def estimate_parameter_errors(parameter_errors, variables, score_vs_variables, threshold, best_variables):
	# Set the maximal acceptable RMSD
	max_score = threshold * score_vs_variables['score'][0]
	# Select the sets of the variables which yeild the RMSD value below the RMSD threshold
	idx = min(range(len(score_vs_variables['score'])), key=lambda i: abs(score_vs_variables['score'][i]-max_score))
	for name in variables:
		var_min = np.amin(score_vs_variables[name][0:idx])
		var_max = np.amax(score_vs_variables[name][0:idx])
		var_opt = best_variables[name]['value']
		var_dev1 = np.absolute(var_opt - var_min)
		var_dev2 = np.absolute(var_opt - var_max)
		var_dev = np.amax([var_dev1, var_dev2])
		parameter_errors[name] = var_dev