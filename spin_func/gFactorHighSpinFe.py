'''
Calculates the g-factor of a high-spin Fe(III)
'''

import numpy as np
from scipy import linalg
from constants import const
from spinMatrices import spinMatrices
from math_func.euler2rotationMatrix import euler2rotationMatrix
import warnings
warnings.filterwarnings("ignore",category=RuntimeWarning)

def gFactorHighSpinFe(spin, mwFreq):
	# Actual g factors
	g_xx = spin['ga'][0]
	g_yy = spin['ga'][1]
	g_zz = spin['ga'][2]
	# g tensor in the molecular frame
	gv = np.array([g_xx, g_yy, g_zz])
	GT_M = np.diag(gv)
	# ZFS constants
	D = spin['D'][0]
	E = spin['D'][1]
	# ZFS tensor in the molecular frame
	Dv = np.array([-D/3 + E, -D/3 - E, 2*D/3])
	DT_M = np.diag(Dv)
	# Electron spin operators
	SZ, SP, SM, SX, SY = spinMatrices(2.5)
	# Calculate effective g factor for 3 orientations along x,y,z
	orient = np.array([[0,         0.5*np.pi, 0],
					   [0.5*np.pi, 0.5*np.pi, 0],
					   [0,                 0, 0]])
	geff = np.zeros(3)
	for i in range(3):
		# Rotation matrix
		RM_L2M = euler2rotationMatrix(orient[i])
		RM_M2L = RM_L2M.T
		# Electron Zeeman term of the Hamiltonian
		GT_L = np.dot(RM_L2M, np.dot(GT_M, RM_M2L))
		GT_L = 0.5 * (GT_L + GT_L.T)
		ez_M = np.array([0.0, 0.0, 1.0])
		ez_L = np.dot(GT_L, ez_M.T)
		Hez = const['Fez'] * ( (ez_L[0] - 1j*ez_L[1]) * 0.5 * SP + \
							   (ez_L[0] + 1j*ez_L[1]) * 0.5 * SM + \
							   ez_L[2] * SZ )
		# ZFS term of the Hamiltonian
		DT_L = np.dot(RM_L2M, np.dot(DT_M, RM_M2L))
		DT_L = 0.5 * (DT_L + DT_L.T)
		Hzfs = DT_L[2][2] * np.dot(SZ, SZ) + \
			   0.25 * (DT_L[0][0] + DT_L[1][1]) * np.dot(SP, SM) + \
			   0.25 * (DT_L[0][0] + DT_L[1][1]) * np.dot(SM, SP) + \
			   0.25 * (DT_L[0][2] - 1j*DT_L[1][2]) * (np.dot(SZ, SP) + np.dot(SP, SZ)) + \
			   0.25 * (DT_L[0][2] + 1j*DT_L[1][2]) * (np.dot(SZ, SM) + np.dot(SM, SZ)) + \
			   0.25 * (DT_L[0][0] - DT_L[1][1] - 1j*(DT_L[0][1] + DT_L[1][0])) * np.dot(SP, SP) + \
			   0.25 * (DT_L[0][0] - DT_L[1][1] + 1j*(DT_L[0][1] + DT_L[1][0])) * np.dot(SM, SM)
		# Calculate resonance fields via the eignenfield method
		E36 = np.eye(36)
		E6 = np.eye(6)
		A = const['GHz2MHz'] * mwFreq * E36 + np.kron(E6, Hzfs.T) - np.kron(Hzfs, E6)
		B = -np.kron(E6, Hez.T) + np.kron(Hez, E6)
		w, vr = linalg.eig(A, B)
		N = 0
		B0v = []
		for k in range(w.size):
			if (np.imag(w[k]) < 0.0001) & (np.real(w[k]) >= 0) & (np.real(w[k]) < 10):
				N += 1
				B0v.append(np.real(w[k]))
		# Find the resonance field of the lowest doublet
		B0 = 0
		dwMin = 0.1
		for k in range(N):
			Htot = B0v[k] * Hez + Hzfs
			w, vr = linalg.eig(Htot)
			w = sorted(np.real(w))
			dw = np.abs(w[0] - w[1]) - const['GHz2MHz'] * mwFreq
			if np.abs(dw) < dwMin:
				dwMin = np.abs(dw)
				B0 = B0v[k]
		# Calculate effective g factor
		geff[i] = const['GHz2MHz'] * mwFreq / (const['Fez'] * B0)
	return geff	