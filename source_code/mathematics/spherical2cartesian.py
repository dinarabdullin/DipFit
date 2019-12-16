'''
Convert spherical coordinates into cartestian coordinates
'''

import numpy as np


def spherical2cartesian(vs):
	N = vs.shape[0]
	vc = np.zeros((N,3))
	for i in range(N):
		r = vs[i][0]
		theta = vs[i][1]
		phi = vs[i][2]
		vc[i][0] = r * np.sin(theta) * np.cos(phi)
		vc[i][1] = r * np.sin(theta) * np.sin(phi)
		vc[i][2] = r * np.cos(theta)
	return vc