'''
Plots of the results of the simulation
'''

import sys
import numpy as np
from simulation.graphics.plot_spectrum import plot_spectrum
from simulation.graphics.plot_spectrum_vs_theta import plot_spectrum_vs_theta
from simulation.graphics.plot_spectrum_vs_par_3d import plot_spectrum_vs_par_3d
from simulation.graphics.plot_spectrum_vs_par_2d import plot_spectrum_vs_par_2d
from simulation.graphics.plot_timetrace import plot_timetrace
from simulation.graphics.plot_pb_vs_temp import plot_pb_vs_temp
from supplement.constants import const

def plot_simulation_data(sim, exp, settings, calc, output):
    sys.stdout.write('Plotting the simulation results... ')
    
    # Plot the dipolar spectrum
    if settings['spc']:
        filename = output['directory'] + 'spc.png'
        plot_spectrum(sim.f, sim.spc, exp['f'], exp['spc'], calc, settings['faxis_normalized'], sim.fn, output['save_figures'], filename)
    
    # Plot the time trace
    if settings['timetrace']:
        filename = output['directory'] + 'timetrace.png'
        plot_timetrace(sim.t, sim.sig, exp['t'], exp['sig'], output['save_figures'], filename)  
    
    # Plot the dipolar spectrum vs theta
    if settings['spc_vs_theta']:
        filename = output['directory'] + 'spc_vs_theta.png'
        plot_spectrum_vs_theta(sim.f, sim.spc, sim.theta_bins, sim.spc_vs_theta, settings['faxis_normalized'], sim.fn, output['save_figures'], filename)
    
    # Plot the dipolar spectrum vs xi
    if settings['spc_vs_xi']:
        filename = output['directory'] + 'spc_vs_xi.png'
        par_label = r'$\mathit{\xi}$ ' + u'(\N{DEGREE SIGN})'
        if (settings['plot_3d']):
            plot_spectrum_vs_par_3d(sim.f, sim.xi_bins, sim.spc_vs_xi, settings['faxis_normalized'], sim.fn, output['save_figures'], filename, par_label, True)
        else:
            plot_spectrum_vs_par_2d(sim.f, sim.xi_bins, sim.spc_vs_xi, settings['faxis_normalized'], sim.fn, output['save_figures'], filename, par_label, False)
    
    # Plot the dipolar spectrum vs phi
    if settings['spc_vs_phi']:
        filename = output['directory'] + 'spc_vs_phi.png'
        par_label = r'$\mathit{\phi}$ '  + u'(\N{DEGREE SIGN})'
        if (settings['plot_3d']):
            plot_spectrum_vs_par_3d(sim.f, sim.phi_bins, sim.spc_vs_phi, settings['faxis_normalized'], sim.fn, output['save_figures'], filename, par_label, False)
        else:
            plot_spectrum_vs_par_2d(sim.f, sim.phi_bins, sim.spc_vs_phi, settings['faxis_normalized'], sim.fn, output['save_figures'], filename, par_label, False)
    
    # Plot the dipolar spectrum vs E/D
    if settings['spc_vs_E']:
        filename = output['directory'] + 'spc_vs_E.png'
        par_label = r'$\mathit{E / D}$'
        if (settings['plot_3d']):
            plot_spectrum_vs_par_3d(sim.f, sim.E_bins, sim.spc_vs_E, settings['faxis_normalized'], sim.fn, output['save_figures'], filename, par_label, False)
        else:
            plot_spectrum_vs_par_2d(sim.f, sim.E_bins, sim.spc_vs_E, settings['faxis_normalized'], sim.fn, output['save_figures'], filename, par_label, False)      

    # Plot the dipolar spectrum vs temperature
    if settings['spc_vs_temp']:    
        filename = output['directory'] + 'spc_vs_temp.png'
        par_label = r'$\mathit{T}$ $(K)$'
        plot_spectrum_vs_par_2d(sim.f, sim.temp_bins, sim.spc_vs_temp, settings['faxis_normalized'], sim.fn, output['save_figures'], filename, par_label, False)
        filename = output['directory'] + 'pb_vs_temp.png'
        plot_pb_vs_temp(sim.g, sim.temp_bins, sim.pb_vs_temp, output['save_figures'], filename)    
    
    sys.stdout.write('[DONE]\n\n')