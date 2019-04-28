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
	
def score_function(var, indices, fixed, sim, exp, spinA, spinB, calc):
    # Set variables
    var_all = set_variables(var, indices, fixed)
    # Calculate dipolar frequencies
    fdd, theta = sim.dipolar_frequencies(var_all, sim.gA, sim.gB, sim.qA, sim.qB, False)
    # Calculate weights 
    if not (indices['temp'] == -1):
        weightB = flip_probabilities(spinB['g'], sim.gB, exp['magnField'], var_all['temp'])
    else:
        weightB = sim.weightB
    # Calculate the dipolar spectrum or the dipolar timetrace and score them against the experimental ones
    if calc['fitted_data'] == 'spectrum':
        # Calculate the dipolar spectrum
        spcSim = sim.dipolar_spectrum(fdd, weightB)
        # Calculate the score
        score = rmsd(spcSim, exp['spc'], exp['f'], calc['f_min'], calc['f_max'])         
    elif calc['fitted_data'] == 'timetrace':
        # Calculate the dipolar timetrace
        #sigSim = sim.dipolar_timetrace(fdd, weightB)
        sigSim = sim.dipolar_timetrace_fast(fdd, weightB)
        # Calculate the score
        score = rmsd(sigSim, exp['sig'], exp['t'], calc['t_min'], calc['t_max'])
    return score
	
def get_fit(var, indices, fixed, sim, exp, spinA, spinB, calc):
    # Set variables
    var_all = set_variables(var, indices, fixed)
    # Calculate dipolar frequencies
    fdd, theta = sim.dipolar_frequencies(var_all, sim.gA, sim.gB, sim.qA, sim.qB, False)
    # Calculate weights 
    if not (indices['temp'] == -1):
        weightB = flip_probabilities(spinB['g'], sim.gB, exp['magnField'], var_all['temp'])
    else:
        weightB = sim.weightB
    # Calculate the dipolar spectrum or the dipolar timetrace
    if calc['fitted_data'] == 'spectrum':
        # Calculate the dipolar spectrum
        fit = sim.dipolar_spectrum(fdd, weightB)        
    elif calc['fitted_data'] == 'timetrace':
        # Calculate the dipolar timetrace
        #fit = sim.dipolar_timetrace(fdd, weightB)
        fit = sim.dipolar_timetrace_fast(fdd, weightB)
    return fit