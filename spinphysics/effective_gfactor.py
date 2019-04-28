'''
Calculates an effective g-factor for a given orientation of magnetic field
'''

import numpy as np

def effective_gfactor(g, eB0):
	# Number of orientations of magnetic field
	N = eB0.shape[0]
	# Calculate effective g factors
	geff = np.zeros(N)
	for i in range(N):
		geff[i] = np.sqrt( (g[0] * eB0[i][0])**2 + 
						   (g[1] * eB0[i][1])**2 + 
						   (g[2] * eB0[i][2])**2 )
	return geff