'''
All valuable constants
'''

import numpy as np


const = {}
const['Hz2MHz'] = 1e-6
const['MHz2Hz'] = 1e6
const['GHz2MHz'] = 1e3
const['mT2T'] = 1e-3
const['T2mT'] = 1e3
const['nm2m'] = 1e-9
const['deg2rad'] = np.pi / 180.0
const['rad2deg'] = 180.0 / np.pi
const['wn2MHz'] = 29979.0
const['plankConstant'] = 6.626070040e-34 # J*s
const['bohrMagneton'] = 9.274009994e-24 # J/T
const['bolzmannConstant'] = 1.38064852e-23 # J/K
const['ge'] = 2.0023 # free electron g factor
const['vacuumPermeabilityOver4Pi'] = 1e-7 # T*m/A
const['Fez'] = const['Hz2MHz'] * const['bohrMagneton'] / const['plankConstant'] # MHz/T
const['Fdd'] = const['Hz2MHz'] * const['vacuumPermeabilityOver4Pi'] * const['bohrMagneton']**2 / (const['plankConstant'] * const['nm2m']**3) # MHz

const['variableNames'] = [
	'r_mean',
	'r_width', 
	'xi_mean', 
	'xi_width', 
	'phi_mean', 
	'phi_width', 
	'temp']
	
const['longVariableNames'] = {
	'r_mean'   : 'r mean (nm)',
	'r_width'  : 'r width (nm)', 
	'xi_mean'  : 'xi mean (deg)',
	'xi_width' : 'xi width (deg)', 
	'phi_mean' : 'phi mean (deg)', 
	'phi_width': 'phi mean (deg)',
	'temp'     : 'temperature (K)'}

const['variableScales'] = {
	'r_mean'   : 1.0,
	'r_width'  : 1.0, 
	'xi_mean'  : const['deg2rad'], 
	'xi_width' : const['deg2rad'], 
	'phi_mean' : const['deg2rad'], 
	'phi_width': const['deg2rad'], 
	'temp'     : 1.0}
	
const['variableLabels'] = {
	'r_mean'   : r'$\langle\mathit{r}\rangle$ (nm)',
	'r_width'  : r'$\mathit{\Delta r}$ (nm)', 
	'xi_mean'  : r'$\langle\mathit{\xi}\rangle$ $^\circ$', 
	'xi_width' : r'$\mathit{\Delta\xi}$ $^\circ$', 
	'phi_mean' : r'$\langle\mathit{\varphi}\rangle$ $^\circ$', 
	'phi_width': r'$\mathit{\Delta\varphi}$ $^\circ$', 
	'temp'     : r'Temperature (K)'}