'''
Score function
'''

import numpy as np
from mathematics.rmsd import rmsd
from spinphysics.flip_probabilities import flip_probabilities
from supplement.constants import const


def set_variables(var, indices, fixed):
	var_all = {}
	for name in const['variableNames']:
		index = indices[name]
		if not (index == -1):
			var_all[name] = var[index]
		else:
			var_all[name] = fixed[name]
	return var_all

	
def score_function(var, fitSettings, simulator, expData, spinA, spinB, calcSettings):
    # Set variables
    indices = fitSettings['variables']['indices']
    fixed = fitSettings['variables']['fixed']
    var_all = set_variables(var, indices, fixed)
    # Calculate dipolar frequencies
    fdd, theta = simulator.dipolar_frequencies(var_all, spinA, spinB)
    # Calculate weights for different g-values of spin B
    pB = []
    if not (indices['temp'] == -1) and calcSettings['g_selectivity'] and (spinB['type'] == "anisotropic"):
        pB = flip_probabilities(simulator.gB, spinB['g'], calcSettings['magnetic_field'], var_all['temp'])
    # Simulate the spectrum or the time trace and score them
    if fitSettings['settings']['fitted_data'] == 'spectrum':
        # Calculate the dipolar spectrum
        spcSim = simulator.dipolar_spectrum(fdd, pB)
        # Calculate the score
        score = rmsd(spcSim, expData['spc'], expData['f'], calcSettings['f_min'], calcSettings['f_max'])         
    elif fitSettings['settings']['fitted_data'] == 'timetrace':
        # Calculate the dipolar timetrace
        #sigSim = simulator.dipolar_timetrace(fdd, pB)
        sigSim = simulator.dipolar_timetrace_fast(fdd, pB)
        # Calculate the score
        score = rmsd(sigSim, expData['sig'], expData['t'], calcSettings['t_min'], calcSettings['t_max'])
    return score

	
def get_fit(var, fitSettings, simulator, expData, spinA, spinB, calcSettings):
    # Set variables
    indices = fitSettings['variables']['indices']
    fixed = fitSettings['variables']['fixed']
    var_all = set_variables(var, indices, fixed)
    # Calculate dipolar frequencies
    fdd, theta = simulator.dipolar_frequencies(var_all, spinA, spinB)
    # Calculate weights for different g-values of spin B
    pB = []
    if not (indices['temp'] == -1) and calcSettings['g_selectivity'] and (spinB['type'] == "anisotropic"):
        pB = flip_probabilities(simulator.gB, spinB['g'], calcSettings['magnetic_field'], var_all['temp'])
    # Calculate the dipolar spectrum or the dipolar time trace
    if fitSettings['settings']['fitted_data'] == 'spectrum':
        # Calculate the dipolar spectrum
        fit = simulator.dipolar_spectrum(fdd, pB)        
    elif fitSettings['settings']['fitted_data'] == 'timetrace':
        # Calculate the dipolar timetrace
        #fit = simulator.dipolar_timetrace(fdd, pB)
        fit = simulator.dipolar_timetrace_fast(fdd, pB)
    return fit