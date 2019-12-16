'''
Estimate the error of the optimized fitting parameters
'''

import sys
import numpy as np
from fitting.generation import Generation
from supplement.constants import const


def calculation_error(best_parameters, fitSettings, simulator, expData, spinA, spinB, calcSettings):
    # Create a generation
    N = 100
    generation = Generation(N)
    generation.first_generation(fitSettings['variables']['bounds'])
    # Set all genes to the optimized values
    for name in const['variableNames']:
        index = fitSettings['variables']['indices'][name]
        if not (index == -1):
            for j in range(N):
                generation.chromosomes[j].genes[index] = best_parameters[name]['value']
    # Score the generation
    generation.score_chromosomes(fitSettings, simulator, expData, spinA, spinB, calcSettings)
    # Sort chromosomes according to their score
    generation.sort_chromosomes()
    # Determine the variation of the score
    score_min = generation.chromosomes[0].score
    score_max = generation.chromosomes[N-1].score
    score_var = score_max - score_min
    return  [score_var, score_min]


def estimate_parameter_errors(parameter_errors, variables, score_vs_variables, score_threshold, best_variables):
    idx = min(range(len(score_vs_variables['score'])), key=lambda i: abs(score_vs_variables['score'][i]-score_threshold))
    for name in variables:
        var_min = np.amin(score_vs_variables[name][0:idx])
        var_max = np.amax(score_vs_variables[name][0:idx])
        var_opt = best_variables[name]['value']
        var_dev1 = np.absolute(var_opt - var_min)
        var_dev2 = np.absolute(var_opt - var_max)
        var_dev = np.amax([var_dev1, var_dev2])
        parameter_errors[name] = var_dev


def calculate_score_vs_parameters(best_parameters, valSettings, fitSettings, simulator, expData, spinA, spinB, calcSettings, score_threshold):
	# Initialize arrays
	score_vs_parameters = []
	parameter_errors = {}
	# Set the number of parameter sets
	M = len(valSettings['variables'])
	# Set the number of Monte-Carlo samples in each calculation
	N = valSettings['Ns']
	# Increment over different sets of variables
	for i in range(M):
		# Number of variables 
		dim = len(valSettings['variables'][i])
		# Create a generation
		generation = Generation(N)
		generation.first_generation(fitSettings['variables']['bounds'])
		# Set all genes to the optimized values except the ones given in varError  
		for name in const['variableNames']:
			index = fitSettings['variables']['indices'][name]
			if not (index == -1) and not (name in valSettings['variables'][i]):
				for j in range(N):
					generation.chromosomes[j].genes[index] = best_parameters[name]['value']
		# Score the generation
		generation.score_chromosomes(fitSettings, simulator, expData, spinA, spinB, calcSettings)
		# Sort chromosomes according to their score
		generation.sort_chromosomes()
		# Store the variables and the corresponding score values
		score_vs_variables = {}
		for name in valSettings['variables'][i]:
			score_vs_variables[name] = []
			index = fitSettings['variables']['indices'][name]
			for j in range(N):
				score_vs_variables[name].append(generation.chromosomes[j].genes[index])
		score_vs_variables['score'] = []
		for j in range(N):
			score_vs_variables['score'].append(generation.chromosomes[j].score)
		score_vs_parameters.append(score_vs_variables)
		# Estimate the parameter errors
		estimate_parameter_errors(parameter_errors, valSettings['variables'][i], score_vs_variables, score_threshold, best_parameters)      
	return [score_vs_parameters, parameter_errors]	

