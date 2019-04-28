'''
Calculate random points on a sphere
'''

import math
import numpy as np
from spherical2cartesian import spherical2cartesian

def rand_points_on_sphere(N):
	# Generate N sets of spherical coordinates
	r = np.ones(N)
	theta = np.arccos(2.0 * np.random.random_sample(N) - 1.0)
	phi = 2.0 * np.pi * np.random.random_sample(N)
	vs = np.array([r, theta, phi])
	vs = vs.T
	# Convert spherical coordinates into Cartesian coordinates
	vc = spherical2cartesian(vs)
	return vc