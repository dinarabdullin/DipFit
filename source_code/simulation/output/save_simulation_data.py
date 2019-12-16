'''
Save the results of the simulation
'''

import sys
from simulation.output.save_spectrum import save_spectrum
from simulation.output.save_spectrum_vs_theta import save_spectrum_vs_theta
from simulation.output.save_spectrum_vs_par import save_spectrum_vs_par
from simulation.output.save_timetrace import save_timetrace
from supplement.constants import const


def save_simulation_data(simulator, simSettings, expData, outputSettings):
    if outputSettings['save_data']:
        sys.stdout.write('Saving the simulation results into the directory:\n')
        sys.stdout.write(outputSettings['directory'])
        if simSettings['modes']['spc']:
            filename = outputSettings['directory'] + 'spc.dat'
            save_spectrum(simulator.f, simulator.spc, expData['f'], expData['spc'], filename)
        if simSettings['modes']['timetrace']:
            filename = outputSettings['directory'] + 'timetrace.dat'
            save_timetrace(simulator.t, simulator.sig, expData['t'], expData['sig'], filename)
            filename = outputSettings['directory'] + 'timetrace_sim.dat'
            save_timetrace(simulator.t, simulator.sig, [], [], filename)
        if simSettings['modes']['spc_vs_theta']:
            filename = outputSettings['directory'] + 'spc_vs_theta.dat'
            save_spectrum_vs_theta(simulator.f, simulator.theta_bins, simulator.spc_vs_theta, filename)
        if simSettings['modes']['spc_vs_xi']:
            filename = outputSettings['directory'] + 'spc_vs_xi.dat'
            save_spectrum_vs_par(simulator.f, simulator.xi_bins, simulator.spc_vs_xi, filename)
        if simSettings['modes']['spc_vs_phi']:
            filename = outputSettings['directory'] + 'spc_vs_phi.dat'
            save_spectrum_vs_par(simulator.f, simulator.phi_bins, simulator.spc_vs_phi, filename) 
        if simSettings['modes']['spc_vs_temp']:
            filename = outputSettings['directory'] + 'spc_vs_temp.dat'
            save_spectrum_vs_par(simulator.f, simulator.temp_bins, simulator.spc_vs_temp, filename)
            filename = outputSettings['directory'] + 'depth_vs_temp.dat'
            save_spectrum_vs_par(simulator.g, simulator.temp_bins, simulator.depth_vs_temp, filename)
        sys.stdout.write(' [DONE]\n\n')