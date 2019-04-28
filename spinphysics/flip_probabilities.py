'''
Weights that are assigned to the different components of the dipolar spectrum 
based on the anisotropic polarization of spin B
'''

import math
import numpy as np
from supplement.constants import const

def flip_probabilities(g, geff, B0, T):
	# Number of different g values
	N = geff.shape[0]
	p = np.zeros(N)
	for i in range(N):
		boltzmann = np.exp(-const['bohrMagneton'] * B0 * geff[i] / (const['bolzmannConstant'] * T))
		p[i] = 2 * boltzmann / ((1.0 + boltzmann) * (1.0 + boltzmann))
	return p