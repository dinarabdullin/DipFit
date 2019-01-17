'''
Calculates the probability that a high-spin Fe3+ spin experiences an odd number of flips in the RIDME interval Tmix
'''

import math
import numpy as np
from constants import const

def flipProbabilities(g, geff, B0, T):
	# Number of different g values
	N = geff.shape[0]
	p = np.zeros(N)
	for i in range(N):
		boltzmann = np.exp(-const['bohrMagneton'] * B0 * geff[i] / (const['bolzmannConstant'] * T))
		p[i] = 2 * boltzmann / ((1.0 + boltzmann) * (1.0 + boltzmann))
	return p