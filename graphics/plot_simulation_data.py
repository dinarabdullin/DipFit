'''
Plots of the results of the simulation
'''

import sys
import numpy as np
from plot_spectrum import plot_spectrum
from plot_spectrum_vs_theta import plot_spectrum_vs_theta
from plot_spectrum_vs_par import plot_spectrum_vs_par
from plot_timetrace import plot_timetrace
from constants import const
from matplotlib import rcParams

def plot_simulation_data(sim, exp, settings, calc, output):
	sys.stdout.write('Plotting the simulation results... ')
	rcParams['font.size'] = 14
	# Plot the dipolar spectrum
	if settings['spc']:
		filename = output['directory'] + 'spc.png'
		plot_spectrum(sim.f, sim.spc, exp['f'], exp['spc'], calc, settings['faxis_normalized'], sim.fn, output['save_figures'], filename)
	# Plot the dipolar spectrum vs theta
	if settings['spc_vs_theta']:
		filename = output['directory'] + 'spc_vs_theta.png'
		plot_spectrum_vs_theta(sim.f, sim.spc, const['theta'], sim.spc_vs_theta, settings['faxis_normalized'], sim.fn, output['save_figures'], filename)
	# Plot the dipolar spectrum vs xi
	if settings['spc_vs_xi']:
		par_label = r'$\mathit{\xi}$ (degree)'
		filename = output['directory'] + 'spc_vs_xi.png'
		plot_spectrum_vs_par(sim.f, const['xi'], sim.spc_vs_xi, settings['faxis_normalized'], sim.fn, output['save_figures'], filename, par_label, True)
	# Plot the dipolar spectrum vs phi
	if settings['spc_vs_phi']:
		par_label = r'$\mathit{\phi}$ (degree)'
		filename = output['directory'] + 'spc_vs_phi.png'
		plot_spectrum_vs_par(sim.f, const['phi'], sim.spc_vs_phi, settings['faxis_normalized'], sim.fn, output['save_figures'], filename, par_label, False)
	# Plot the dipolar spectrum vs E/D
	if settings['spc_vs_E']:
		par_label = r'$\mathit{E / D}$'
		filename = output['directory'] + 'spc_vs_E.png'
		plot_spectrum_vs_par(sim.f, const['E'], sim.spc_vs_E, settings['faxis_normalized'], sim.fn, output['save_figures'], filename, par_label, False)        
	# Plot the time trace
	if settings['timetrace']:
		filename = output['directory'] + 'timetrace.png'
		plot_timetrace(sim.t, sim.sig, exp['t'], exp['sig'], output['save_figures'], filename)    
	sys.stdout.write('[DONE]\n\n')