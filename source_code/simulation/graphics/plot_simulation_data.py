'''
Plot of the results of the simulation
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


def plot_simulation_data(simulator, simSettings, expData, calcSettings, outputSettings):
    sys.stdout.write('Plotting the simulation results... ')
    if simSettings['modes']['spc']:
        filename = outputSettings['directory'] + 'spc.png'
        plot_spectrum(simulator.f, simulator.spc, expData['f'], expData['spc'], calcSettings, 
                      simSettings['settings']['faxis_normalized'], simulator.fn, 
                      outputSettings['save_figures'], filename)
    if simSettings['modes']['timetrace']:
        filename = outputSettings['directory'] + 'timetrace.png'
        plot_timetrace(simulator.t, simulator.sig, expData['t'], expData['sig'], 
                       outputSettings['save_figures'], filename)  
    if simSettings['modes']['spc_vs_theta']:
        filename = outputSettings['directory'] + 'spc_vs_theta.png'
        plot_spectrum_vs_theta(simulator.f, simulator.spc, simulator.theta_bins, simulator.spc_vs_theta, calcSettings, 
                               simSettings['settings']['faxis_normalized'], simulator.fn, 
                               outputSettings['save_figures'], filename)
    if simSettings['modes']['spc_vs_xi']:
        filename = outputSettings['directory'] + 'spc_vs_xi.png'
        par_label = r'$\mathit{\xi}$ ' + u'(\N{DEGREE SIGN})'
        if simSettings['settings']['plot_3d']:
            plot_spectrum_vs_par_3d(simulator.f, simulator.xi_bins, simulator.spc_vs_xi, calcSettings,
                                    simSettings['settings']['faxis_normalized'], simulator.fn, 
                                    outputSettings['save_figures'], filename, par_label, True)
        else:
            plot_spectrum_vs_par_2d(simulator.f, simulator.xi_bins, simulator.spc_vs_xi, calcSettings, 
                                    simSettings['settings']['faxis_normalized'], simulator.fn, 
                                    outputSettings['save_figures'], filename, par_label, False)
    if simSettings['modes']['spc_vs_phi']:
        filename = outputSettings['directory'] + 'spc_vs_phi.png'
        par_label = r'$\mathit{\phi}$ '  + u'(\N{DEGREE SIGN})'
        if simSettings['settings']['plot_3d']:
            plot_spectrum_vs_par_3d(simulator.f, simulator.phi_bins, simulator.spc_vs_phi, calcSettings, 
                                    simSettings['settings']['faxis_normalized'], simulator.fn, 
                                    outputSettings['save_figures'], filename, par_label, False)
        else:
            plot_spectrum_vs_par_2d(simulator.f, simulator.phi_bins, simulator.spc_vs_phi, calcSettings, 
                                    simSettings['settings']['faxis_normalized'], simulator.fn, 
                                    outputSettings['save_figures'], filename, par_label, False)   
    if simSettings['modes']['spc_vs_temp']:
        filename = outputSettings['directory'] + 'spc_vs_temp.png'
        par_label = r'$\mathit{Temperature (K)}$'
        if simSettings['settings']['plot_3d']:
            plot_spectrum_vs_par_3d(simulator.f, simulator.temp_bins, simulator.spc_vs_temp, calcSettings, 
                                    simSettings['settings']['faxis_normalized'], simulator.fn, 
                                    outputSettings['save_figures'], filename, par_label, False)
        else:
            plot_spectrum_vs_par_2d(simulator.f, simulator.temp_bins, simulator.spc_vs_temp, calcSettings, 
                                    simSettings['settings']['faxis_normalized'], simulator.fn, 
                                    outputSettings['save_figures'], filename, par_label, False)  
    sys.stdout.write('[DONE]\n\n')