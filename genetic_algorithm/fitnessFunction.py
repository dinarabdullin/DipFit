'''
Fitness function
'''

import numpy as np
from copy import deepcopy
from constants import const
from math_func.rmsd import rmsd
from spin_func.flipProbabilities import flipProbabilities

def set_variables(var, indices, fixed):
	var_all = {}
	for name in const['variableNames']:
		index = indices[name]
		if not (index == -1):
			var_all[name] = var[index]
		else:
			var_all[name] = fixed[name]
	return var_all
	
def fitness_function(var, indices, fixed, sim, exp, spinA, spinB, calc):
    # Set variables
    var_all = set_variables(var, indices, fixed)
    # Calculate dipolar frequencies
    fdd, theta = sim.dipolar_frequencies(var_all, sim.gA, sim.gB, sim.qA, sim.qB, False)
    # Calculate the dipolar spectrum
    if not (indices['temp'] == -1):
        # Calculate weights 
        weightB = flipProbabilities(spinB['g'], sim.gB, exp['magnField'], var_all['temp'])
        spcSim = sim.dipolar_spectrum(fdd, weightB)
    else:
        spcSim = sim.dipolar_spectrum(fdd)
    # Calculate RMSD between the simulated and experimental spectra
    fitness = rmsd(spcSim, exp['spc'], exp['f'], calc['f_min'], calc['f_max'])
    return fitness
	
def get_fit(var, indices, fixed, sim, exp, spinA, spinB):
    # Set variables
    var_all = set_variables(var, indices, fixed)
    # Calculate dipolar frequencies
    fdd, theta = sim.dipolar_frequencies(var_all, sim.gA, sim.gB, sim.qA, sim.qB, False)
    # Calculate the dipolar spectrum
    if not (indices['temp'] == -1):
        # Calculate weights 
        weightB = flipProbabilities(spinB['g'], sim.gB, exp['magnField'], var_all['temp'])
        spcSim = sim.dipolar_spectrum(fdd, weightB)
    else:
        spcSim = sim.dipolar_spectrum(fdd)
    return spcSim