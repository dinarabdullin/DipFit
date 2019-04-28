'''
Saves the results of the simulation
'''

import sys
from simulation.output.save_spectrum import save_spectrum
from simulation.output.save_spectrum_vs_theta import save_spectrum_vs_theta
from simulation.output.save_spectrum_vs_par import save_spectrum_vs_par
from simulation.output.save_timetrace import save_timetrace
from supplement.constants import const

def save_simulation_data(sim, exp, settings, output):
	if output['save_data']:
		# Output status
		sys.stdout.write('Saving the simulation results into the directory:\n')
		sys.stdout.write(output['directory'])
        
		# Save the dipolar spectrum
		if settings['spc']:
			filename = output['directory'] + 'spc.dat'
			save_spectrum(sim.f, sim.spc, exp['f'], exp['spc'], filename)
		
        # Save the dipolar spectrum vs theta
		if settings['spc_vs_theta']:
			filename = output['directory'] + 'spc_vs_theta.dat'
			save_spectrum_vs_theta(sim.f, sim.theta_bins, sim.spc_vs_theta, filename)
		
        # Save the dipolar spectrum vs xi
		if settings['spc_vs_xi']:
			filename = output['directory'] + 'spc_vs_xi.dat'
			save_spectrum_vs_par(sim.f, sim.xi_bins, sim.spc_vs_xi, filename)
		
        # Save the dipolar spectrum vs phi
		if settings['spc_vs_phi']:
			filename = output['directory'] + 'spc_vs_phi.dat'
			save_spectrum_vs_par(sim.f, sim.phi_bins, sim.spc_vs_phi, filename)
		
        # Save the dipolar spectrum vs E/D
		if settings['spc_vs_E']:
			filename = output['directory'] + 'spc_vs_E.dat'
			save_spectrum_vs_par(sim.f, sim.E_bins, sim.spc_vs_E, filename)        
		
        # Save the time trace
		if settings['timetrace']:
			filename = output['directory'] + 'timetrace.dat'
			save_timetrace(sim.t, sim.sig, exp['t'], exp['sig'], filename)
			filename = output['directory'] + 'timetrace_sim.dat'
			save_timetrace(sim.t, sim.sig, [], [], filename)
		
        # Output status
		sys.stdout.write(' [DONE]\n\n')