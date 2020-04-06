'''
Estimate the error of the optimized fitting parameters
'''

import sys
import numpy as np
from scipy.optimize import curve_fit
from functools import partial
from fitting.generation import Generation
from supplement.constants import const


def calculate_numerical_error(best_parameters, err_settings, fit_settings, simulator, exp_data, spinA, spinB, calc_settings):
    # Create a generation
    Ns = err_settings['Ns']
    generation = Generation(Ns)
    generation.first_generation(fit_settings['parameters']['bounds'])
    # Set all genes to the optimized values
    for name in const['variable_names']:
        index = fit_settings['parameters']['indices'][name]
        if not (index == -1):
            for j in range(Ns):
                generation.chromosomes[j].genes[index] = best_parameters[name]['value']
    # Score the generation
    generation.score_chromosomes(fit_settings, simulator, exp_data, spinA, spinB, calc_settings)
    # Sort chromosomes according to their score
    generation.sort_chromosomes()
    # Determine the variation of the score
    score_min = generation.chromosomes[0].score
    score_max = generation.chromosomes[Ns-1].score
    numerical_error = score_max - score_min
    sys.stdout.write("Numerical error = %f\n" % (numerical_error))
    sys.stdout.write("Minimal RMSD = %f\n" % (score_min))
    return [numerical_error, score_min]


def calculate_score_threshold(best_score, threshold, numerical_error):
    score_threshold = best_score * threshold
    if (numerical_error > (score_threshold - best_score)):
        score_threshold = best_score + numerical_error
        sys.stdout.write('Warning: The RMSD threshold is too low!\n')
        sys.stdout.write('The RMSD threshold is set to the sum of the minimal RMSD and the numerical error.\n')
    sys.stdout.write("RMSD threshold = %f\n" % (score_threshold))
    return score_threshold


def calculate_score_vs_parameters(best_parameters, err_settings, fit_settings, simulator, exp_data, spinA, spinB, calc_settings):
    sys.stdout.write("Calculating the RMSD in dependence of fitting parameters ...\n")
    score_vs_parameters = []
    Ne = len(err_settings['variables'])
    Ns = err_settings['Ns']
    for i in range(Ne):
        sys.stdout.write('\r')
        sys.stdout.write("Calculation step %d / %d" % (i+1, Ne))
        sys.stdout.flush()
        # Create a generation
        generation = Generation(Ns)
        generation.first_generation(fit_settings['parameters']['bounds'])
        # Set all genes to the optimized values except the ones given in varError  
        for name in const['variable_names']:
            index = fit_settings['parameters']['indices'][name]
            if not (index == -1) and not (name in err_settings['variables'][i]):
                for j in range(Ns):
                    generation.chromosomes[j].genes[index] = best_parameters[name]['value']
        # Score the generation
        generation.score_chromosomes(fit_settings, simulator, exp_data, spinA, spinB, calc_settings)
        # Sort chromosomes according to their score
        generation.sort_chromosomes()
        # Store the variables and the corresponding score values
        score_vs_variables = {}
        for name in err_settings['variables'][i]:
            score_vs_variables[name] = []
            index = fit_settings['parameters']['indices'][name]
            for j in range(Ns):
                score_vs_variables[name].append(generation.chromosomes[j].genes[index])
        score_vs_variables['score'] = []
        for j in range(Ns):
            score_vs_variables['score'].append(generation.chromosomes[j].score)
        score_vs_parameters.append(score_vs_variables)
    sys.stdout.write('\n')
    return score_vs_parameters


def calculate_parameter_errors(variables, score_vs_parameters, best_parameters, score_threshold):
    sys.stdout.write("Calculating the errors of fitting parameters ...\n")
    parameter_errors = {}
    Ne = len(score_vs_parameters)
    for i in range(Ne):
        for name in variables[i]:
            variable_values = np.array(score_vs_parameters[i][name])
            score_values = np.array(score_vs_parameters[i]['score'])
            best_parameter = best_parameters[name]['value']
            parameter_error = calculate_parameter_error(variable_values, score_values, best_parameter, score_threshold)
            if name in parameter_errors:
                if not np.isnan(parameter_error) and not np.isnan(parameter_errors[name]):
                    if (parameter_error > parameter_errors[name]):
                        parameter_errors[name] = parameter_error
            else:
                parameter_errors[name] = parameter_error
    return parameter_errors


def calculate_parameter_error(x_data, y_data, x_opt, threshold):
    Ns = x_data.size
    # Determine the minimal and maximal values of x
    x_min = np.amin(x_data)
    x_max = np.amax(x_data)
    # Set the optimal values of x and y
    if np.isnan(x_opt):
        y_opt = np.amin(y_data)
        idx_y_opt = np.argmin(y_data)
        x_opt = x_data[idx_y_opt]
    else:
        idx_x_opt = min(range(len(x_data)), key=lambda i: abs(x_data[i]-x_opt))
        y_opt = y_data[idx_x_opt]
    # Sort x in ascending order
    x_sorted, y_sorted = zip(*sorted(zip(x_data, y_data)))
    # Determine the uncertainty ranges of x
    idx_x_selected = []
    for i in range(Ns):
        if y_sorted[i] < threshold:
            idx_x_selected.append(i)
    x_left = x_sorted[idx_x_selected[0]]
    x_right = x_sorted[idx_x_selected[-1]]
    # Determine the error of x_opt
    x_dev_left = abs(x_opt - x_left)
    x_dev_right = abs(x_opt - x_right)
    x_error = np.amax([x_dev_left, x_dev_right])
    if (x_left == x_min) and (x_right == x_max):
        x_error = np.nan
    return x_error