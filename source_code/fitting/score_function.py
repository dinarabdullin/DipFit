'''
Score function
'''

import numpy as np
from mathematics.rmsd import rmsd
from spinphysics.flip_probabilities import flip_probabilities
from supplement.constants import const


def set_parameters(parameters, indices, fixed):
	var_all = {}
	for name in const['variable_names']:
		index = indices[name]
		if not (index == -1):
			var_all[name] = parameters[index]
		else:
			var_all[name] = fixed[name]
	return var_all

	
def score_function(parameters, fit_settings, simulator, exp_data, spinA, spinB, calc_settings):
    # Set parameters
    indices = fit_settings['parameters']['indices']
    fixed = fit_settings['parameters']['fixed']
    all_parameters = set_parameters(parameters, indices, fixed)
    # Calculate dipolar frequencies
    fdd, theta = simulator.dipolar_frequencies(all_parameters, spinA, spinB)
    # Calculate weights for different g-values of spin B
    pB = []
    if (not (indices['temp'] == -1) and calc_settings['g_selectivity'] and (spinB['type'] == "anisotropic")):
        pB = flip_probabilities(simulator.gB, spinB['g'], calc_settings['magnetic_field'], all_parameters['temp'])
    # Simulate the spectrum or the time trace and score them
    if fit_settings['settings']['fitted_data'] == 'spectrum':
        # Calculate the dipolar spectrum
        fit = simulator.dipolar_spectrum(fdd, pB)
        # Calculate the score
        score = rmsd(fit, exp_data['spc'], exp_data['f'], calc_settings['f_min'], calc_settings['f_max'])         
    elif fit_settings['settings']['fitted_data'] == 'timetrace':
        # Calculate the dipolar timetrace
        #fit = simulator.dipolar_timetrace(fdd, pB)
        fit = simulator.dipolar_timetrace_fast(fdd, pB)
        # Calculate the score
        score = rmsd(fit, exp_data['sig'], exp_data['t'], calc_settings['t_min'], calc_settings['t_max'])
    return score

	
def get_fit(parameters, fit_settings, simulator, exp_data, spinA, spinB, calc_settings):
    # Set parameters
    indices = fit_settings['parameters']['indices']
    fixed = fit_settings['parameters']['fixed']
    all_parameters = set_parameters(parameters, indices, fixed)
    # Calculate dipolar frequencies
    fdd, theta = simulator.dipolar_frequencies(all_parameters, spinA, spinB)
    # Calculate weights for different g-values of spin B
    pB = []
    if (not (indices['temp'] == -1) and calc_settings['g_selectivity'] and (spinB['type'] == "anisotropic")):
        pB = flip_probabilities(simulator.gB, spinB['g'], calc_settings['magnetic_field'], all_parameters['temp'])
    # Calculate the dipolar spectrum or the dipolar time trace
    if fit_settings['settings']['fitted_data'] == 'spectrum':
        # Calculate the dipolar spectrum
        fit = simulator.dipolar_spectrum(fdd, pB)        
    elif fit_settings['settings']['fitted_data'] == 'timetrace':
        # Calculate the dipolar timetrace
        #fit = simulator.dipolar_timetrace(fdd, pB)
        fit = simulator.dipolar_timetrace_fast(fdd, pB)
    return fit