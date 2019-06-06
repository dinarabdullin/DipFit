'''
Calculate the quantisation axis of a spin center for a given orientation of magnetic field
'''

import numpy as np


def quantisation_axis(g, geff, eB0):
	# Number of different orientation sof magnetic field
	N = eB0.shape[0]
	# Initialize the array with cartesian coordinates
	qa = np.zeros((N,3))
	for i in range(N):
		qa[i][0] = (g[0] / geff[i])**2 * eB0[i][0]
		qa[i][1] = (g[1] / geff[i])**2 * eB0[i][1]
		qa[i][2] = (g[2] / geff[i])**2 * eB0[i][2]
	return qa