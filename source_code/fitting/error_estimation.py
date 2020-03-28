'''
Estimate the error of the optimized fitting parameters
'''

import sys
import numpy as np
from scipy.optimize import curve_fit
from functools import partial
from fitting.generation import Generation
from supplement.constants import const



def calculate_numerical_error(best_parameters, valSettings, fitSettings, simulator, expData, spinA, spinB, calcSettings):
    # Create a generation
    N = valSettings['Ns']
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
    return  score_var


def calculate_score_vs_parameters(best_parameters, valSettings, fitSettings, simulator, expData, spinA, spinB, calcSettings):
	score_vs_parameters = []
	M = len(valSettings['variables'])
	N = valSettings['Ns']
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
	return score_vs_parameters


def calculate_parameter_errors(score_vs_variables, best_variables, valSettings, numerical_error):
    parameter_errors = {}
    M = len(valSettings['variables'])
    for i in range(M):
        for name in valSettings['variables'][i]:
            var = np.array(score_vs_variables[i][name])
            score = np.array(score_vs_variables[i]['score'])
            threshold = valSettings['threshold']
            var_opt = best_variables[name]['value']
            parameter_error = calculate_parameter_error(var, score, name, threshold, var_opt, numerical_error)
            if name in parameter_errors:
                if not np.isnan(parameter_error) and not np.isnan(parameter_errors[name]):
                    if (parameter_error > parameter_errors[name]):
                        parameter_errors[name] = parameter_error
            else:
                parameter_errors[name] = parameter_error
    return parameter_errors


def gauss_fit(x, a, dx, x0, y0):
    return y0 + a * (1 - np.exp(-0.5 * ((x-x0) / dx)**2))


def calculate_parameter_error(x_data, y_data, x_name, threshold, x_opt, numerical_error):
    # Determine the minimal and maximal values of x
    x_min = np.amin(x_data)
    x_max = np.amax(x_data)
    # Create the new x-axis
    Nx = 101
    x_data_interpolated = np.linspace(x_min, x_max, Nx)
    x_data_inc = x_data_interpolated[1] - x_data_interpolated[0]
    x_data_interpolated_lower = x_data_interpolated - x_data_inc * np.ones(x_data_interpolated.size)
    x_data_interpolated_lower[0] = x_data_interpolated[0]
    x_data_interpolated_upper = x_data_interpolated + x_data_inc * np.ones(x_data_interpolated.size)
    x_data_interpolated_upper[-1] = x_data_interpolated[-1]
    # Select the minimal value of y for each x value
    y_data_interpolated = np.zeros(x_data_interpolated.size)
    for i in range(Nx):
        y_data_interpolated[i] = np.amin(y_data[(x_data_interpolated_lower[i] < x_data) & (x_data < x_data_interpolated_upper[i])])  
    # Set the optimal values of x and y
    if np.isnan(x_opt):
        y_opt = np.amin(y_data_interpolated)
        idx_y_opt = np.argmin(y_data_interpolated)
        x_opt = x_data_interpolated[idx_y_opt]
    else:
        idx_x_opt = min(range(len(x_data_interpolated)), key=lambda i: abs(x_data_interpolated[i]-x_opt))
        y_opt = y_data_interpolated[idx_x_opt]
    # Fit the dependence y(x) using the parameters a and x0 of the function gauss_fit:
    func = partial(gauss_fit, x0=x_opt, y0=y_opt)
    popt, pcov = curve_fit(func, x_data_interpolated, y_data_interpolated)
    # Determine the confidence interval
    dx_opt = abs(popt[1])
    x_error = threshold * dx_opt
    x_left = x_opt - x_error
    x_right = x_opt + x_error
    if a <= 0.2 * np.amin(y_data_interpolated):
        x_left = x_min
        x_right = x_max
        x_error_final = np.nan
        y_dev = np.nan
    else:
        if (x_left >= x_min) and (x_right <= x_max):
            x_error_final = x_error
            idx_left = min(range(len(x_data_interpolated)), key=lambda i: abs(x_data_interpolated[i]-x_left))
            idx_right = min(range(len(x_data_interpolated)), key=lambda i: abs(x_data_interpolated[i]-x_right))
            y_left = y_data_interpolated[idx_left]
            y_right = y_data_interpolated[idx_right]
            y_dev_left = abs(y_opt - y_left)
            y_dev_right = abs(y_opt - y_right)
            y_dev = np.amin([y_dev_left, y_dev_right])
        elif (x_left < x_min) and (x_right <= x_max):
            x_left = x_min
            x_error_final = x_error
            idx_right = min(range(len(x_data_interpolated)), key=lambda i: abs(x_data_interpolated[i]-x_right))
            y_right = y_data_interpolated[idx_right]
            y_dev = abs(y_opt - y_right)
        elif (x_left >= x_min) and (x_right > x_max):
            x_right = x_max
            x_error_final = x_error
            idx_left = min(range(len(x_data_interpolated)), key=lambda i: abs(x_data_interpolated[i]-x_left))
            y_left = y_data_interpolated[idx_left]
            y_dev = abs(y_opt - y_left)
        elif (x_left < x_min) and (x_right > x_max):
            x_left = x_min
            x_right = x_max
            x_error_final = np.nan
            y_dev = np.nan
    # Check that the numerical error is inside the RMSD threshold
    if not np.isnan(numerical_error) and not np.isnan(y_dev):
        if numerical_error > y_dev:
            sys.stdout.write('\nWarning: The RMSD threshold is under the numerical error! Increase the threshold.')
    return x_error_final