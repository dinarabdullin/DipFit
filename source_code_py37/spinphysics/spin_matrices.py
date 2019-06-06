'''
Calculate the spin matrices for a given S value
'''

import numpy as np


def spin_matrices(s):
	# Multiplicity
	M = int(2*s + 1)
	# Initialize the spin matrices
	SZ = np.zeros((M,M))
	SP = np.zeros((M,M))
	SM = np.zeros((M,M))
	SX = np.zeros((M,M))
	SY = np.zeros((M,M))
	# Calculate the spin matrices
	SZ = np.diag(np.linspace(-s, s, M)) 
	for i in range(M-1):
		x = np.sqrt(float((i+1) * (M - (i+1))))
		SP[i][i+1] = x
		SM[i+1][i] = x
		SX[i][i+1] = 0.5 * x
		SX[i+1][i] = 0.5 * x
		SY[i][i+1] = -0.5 * (i+1) * x
		SY[i+1][i] = 0.5 * (i+1) * x
	return SZ, SP, SM, SX, SY