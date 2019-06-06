'''
Read the input data from the config file
'''

import io
import sys
import libconf
import numpy as np
from input.load_spectrum import load_spectrum
from input.load_timetrace import load_timetrace
from input.load_solution import load_solution
from input.symmetric_boundaries import symmetric_boundaries
from mathematics.find_max_on_interval import find_max_on_interval
from spinphysics.gfactor_hs_iron import gfactor_hs_iron
from supplement.constants import const


def read_experimental_parameters(config):
    # Init the data container
    exp = {}  
    exp['t'] = []
    exp['sig'] = []
    exp['f'] = []
    exp['spc'] = []
    # Read out the path to the RIDME spectrum
    exp['path_spectrum'] = config.path_spectrum
    # Read out the path to the RIDME time trace
    exp['path_timetrace'] = config.path_timetrace
    if exp['path_spectrum']:
        # Read out the RIDME spectrum
        f, spc = load_spectrum(exp['path_spectrum'])
        exp['f'], exp['spc'] = symmetric_boundaries(f, spc)
    if exp['path_timetrace']:
        # Read out the RIDME time trace
        exp['t'], exp['sig'] = load_timetrace(exp['path_timetrace'])
    # Read out the value of microwave frequency
    exp['mwFreq'] = float(config.mw_frequency)
    # Read out the value of magnetic field
    exp['magnField'] = float(config.magnetic_field)
    return exp


def read_spin_parameters(config, exp):
    # Init the data container
    spinA = {}
    spinB = {}
    # Read out the type of spin center for spin A
    spinA['label'] = config.spinA.label
    # Read out the zero field splitting tensor of spin A
    spinA['D'] = const['wn2MHz'] * np.array([float(config.spinA.D[0]), float(config.spinA.D[1])])
    # Read out the g-tensor of spin A
    spinA['ga'] = np.array([float(config.spinA.g[0]), 
                            float(config.spinA.g[1]), 
                            float(config.spinA.g[2])])
    # Read out the orientation of g-frame relative to the ZFS frame for spin A
    spinA['gFrame'] = const['deg2rad'] * np.array([float(config.spinA.gFrame[0]), 
                                                   float(config.spinA.gFrame[1]), 
                                                   float(config.spinA.gFrame[2])])
    # Calculate the principal values of g-factor for spin A
    if (spinA['label'] == 'hs_fe'):
       raise ValueError('Error: Only spin B can be a high-spin Fe(III). Spin A should be an isotropic spin-1/2!')
       sys.exit(1)
    else:
        spinA['g'] = spinA['ga']
    # Read out the type of spin center for spin B
    spinB['label'] = config.spinB.label
    # Read out the zero field splitting tensor of spin B
    spinB['D'] = const['wn2MHz'] * np.array([float(config.spinB.D[0]), float(config.spinB.D[1])])
    # Read out the g-tensor of spin B
    spinB['ga'] = np.array([float(config.spinB.g[0]), 
                            float(config.spinB.g[1]), 
                            float(config.spinB.g[2])])
    # Read out the orientation of g-frame relative to the ZFS frame for spin B
    spinB['gFrame'] = const['deg2rad'] * np.array([float(config.spinB.gFrame[0]), 
                                                   float(config.spinB.gFrame[1]), 
                                                   float(config.spinB.gFrame[2])])
    # Calculate the principal values of g-factor for spin B
    if (spinB['label'] == 'hs_fe'):
        spinB['g'] = gfactor_hs_iron(spinB, exp['mwFreq'])
    else:
        spinB['g'] = spinB['ga']
    return [spinA, spinB]

    
def read_calculation_mode(config):
    mode = {}
    # Read out the calculation mode
    switch = int(config.mode)
    if (switch == 0):
        mode['simulation'] = 1
        mode['fitting'] = 0
        mode['validation'] = 0
    elif (switch == 1):
        mode['simulation'] = 0
        mode['fitting'] = 1
        mode['validation'] = 0
    elif (switch == 2):
        mode['simulation'] = 0
        mode['fitting'] = 0
        mode['validation'] = 1
    return mode

def read_simulation_parameters(config, exp):
    simPar = {}
    # Simulation settings
    simPar['settings'] = {}
    simPar['settings']['spc'] = int(config.simulation_settings.spc)
    simPar['settings']['spc_vs_theta'] = int(config.simulation_settings.spc_vs_theta) 
    simPar['settings']['spc_vs_xi'] = int(config.simulation_settings.spc_vs_xi)
    simPar['settings']['spc_vs_phi'] = int(config.simulation_settings.spc_vs_phi)
    simPar['settings']['spc_vs_E'] = int(config.simulation_settings.spc_vs_E)
    simPar['settings']['spc_vs_temp'] = int(config.simulation_settings.spc_vs_temp)
    simPar['settings']['timetrace'] = int(config.simulation_settings.timetrace)
    if (exp['f'] == []):
        simPar['settings']['faxis_normalized'] = int(config.simulation_settings.faxis_normalized)
    else:
        simPar['settings']['faxis_normalized'] = 0
    simPar['settings']['plot_3d'] = int(config.simulation_settings.plot_3d)  
    # Simulation variables
    simPar['variables'] = {}
    simPar['variables']['r_mean'] = float(config.simulation_variables.r_mean)
    simPar['variables']['r_width'] = float(config.simulation_variables.r_width)
    simPar['variables']['xi_mean'] = const['deg2rad'] * float(config.simulation_variables.xi_mean)
    simPar['variables']['xi_width'] = const['deg2rad'] * float(config.simulation_variables.xi_width)
    simPar['variables']['phi_mean'] = const['deg2rad'] * float(config.simulation_variables.phi_mean)
    simPar['variables']['phi_width'] = const['deg2rad'] * float(config.simulation_variables.phi_width)
    simPar['variables']['temp'] = float(config.simulation_variables.temp)
    return simPar


def read_fitting_variables(object):
	indices = {}
	bounds = []
	fixed = {}
	count = 0
	for name in const['variableNames']:
		attribute = getattr(object, name)
		vopt = int(attribute.opt)
		if vopt:
			indices[name] = count
			count += 1
			vmin = float(attribute.range[0]) * const['variableScales'][name]
			vmax = float(attribute.range[1]) * const['variableScales'][name]
			bounds.append([vmin, vmax])
		else:
			indices[name] = -1
			fixed[name] = float(attribute.value) * const['variableScales'][name]
	return [indices, bounds, fixed, count]


def read_fitting_parameters(config):
    fitPar = {}
    # Fitting settings
    fitPar['settings'] = {}
    fitPar['settings']['num_generations'] = int(config.genetic_algorithm.num_generations)
    fitPar['settings']['size_generation'] = int(config.genetic_algorithm.size_generation)
    fitPar['settings']['prob_crossover'] = float(config.genetic_algorithm.prob_crossover)
    fitPar['settings']['prob_mutation'] = float(config.genetic_algorithm.prob_mutation)
    fitPar['settings']['display_graphics'] = int(config.genetic_algorithm.display_graphics)
    # Fitting variables
    fitPar['variables'] = {}
    fitPar['variables']['indices'] = {}
    fitPar['variables']['bounds'] = []
    fitPar['variables']['fixed'] = {}
    fitPar['variables']['indices'], fitPar['variables']['bounds'], fitPar['variables']['fixed'], fitPar['variables']['size'] = read_fitting_variables(config.fitting_variables)
    # Error estimates
    fitPar['errors'] = {}
    fitPar['errors']['variables'] = config.error_estimation.variables
    fitPar['errors']['Ns'] = int(config.error_estimation.Ns)
    fitPar['errors']['threshold'] = 0.01 * float(config.error_estimation.threshold)
    return fitPar


def read_calculation_parameters(config, mode, exp, simPar):
    calc = {}
    solution = {}
    # General calculation parameters
    calc['Ns'] = int(config.calculation_settings.Ns) 
    calc['r_distr'] = config.calculation_settings.r_distr
    calc['xi_distr'] = config.calculation_settings.xi_distr
    calc['phi_distr'] = config.calculation_settings.phi_distr
    calc['f_min'] = float(config.calculation_settings.fmin)
    calc['f_max'] = float(config.calculation_settings.fmax)
    calc['t_min'] = float(config.calculation_settings.tmin)
    calc['t_max'] = float(config.calculation_settings.tmax)
    # Set the maximal amplitude of the RIDME spectrum and the maximal dipolar frequency
    if (exp['f'] == []):
        calc['spc_max'] = 1.0
    else:
        if calc['f_max']:
            calc['spc_max'] = find_max_on_interval(exp['spc'],exp['f'],calc['f_min'],calc['f_max'])   
        else:
            calc['f_max'] = np.amax(exp['f'])
            calc['spc_max'] = find_max_on_interval(exp['spc'],exp['f'],calc['f_min'],calc['f_max'])
    # Set the maximal amplitude of the RIDME time trace and the maximal time point    
    if not (exp['t'] == []):
            if not calc['t_max']:
                calc['t_max'] = np.amax(exp['t'])
    # The parameters used in simulation
    if mode['simulation']:
        if simPar['settings']['spc_vs_theta']:
            calc['theta_bins'] = np.linspace(float(config.calculation_settings.theta_val[0]), 
                                             float(config.calculation_settings.theta_val[1]), 
                                             int(config.calculation_settings.theta_val[2]))
        else:
            calc['theta_bins'] = []
        if simPar['settings']['spc_vs_xi']:
            calc['xi_bins'] = np.linspace(float(config.calculation_settings.xi_val[0]), 
                                          float(config.calculation_settings.xi_val[1]), 
                                          int(config.calculation_settings.xi_val[2]))
        else:
            calc['xi_bins'] = []
        if simPar['settings']['spc_vs_phi']:
            calc['phi_bins'] = np.linspace(float(config.calculation_settings.phi_val[0]), 
                                           float(config.calculation_settings.phi_val[1]), 
                                           int(config.calculation_settings.phi_val[2]))
        else:
            calc['phi_bins'] = []
        if simPar['settings']['spc_vs_E']:
            calc['E_bins'] = np.linspace(float(config.calculation_settings.E_val[0]), 
                                         float(config.calculation_settings.E_val[1]), 
                                         int(config.calculation_settings.E_val[2]))
        else:
            calc['E_bins'] = []
        if simPar['settings']['spc_vs_temp']:
            temp_array = config.calculation_settings.temp_val
            temp_bins = []
            for temp in temp_array:
                temp_bins.append(float(temp))    
            calc['temp_bins'] = np.array(temp_bins)
        else:
            calc['temp_bins'] = []
    else:
        calc['theta_bins'] = []
        calc['xi_bins'] = []
        calc['phi_bins'] = []
        calc['E_bins'] = []
        calc['temp_bins'] = []
    # The parameters used in fitting
    if mode['fitting'] or mode['validation']:
        calc['fitted_data'] = config.calculation_settings.fitted_data
        if (calc['fitted_data'] == 'spectrum') and (exp['f'] == []):
            raise ValueError('The dipolar spectrum was not specified!')
            sys.exit(1)
        elif (calc['fitted_data'] == 'timetrace') and (exp['t'] == []):
            raise ValueError('The dipolar time trace was not specified!')
            sys.exit(1) 
    # The parameters used in validation
    if mode['validation']:
        path_solution = config.calculation_settings.path_solution
        solution = load_solution(path_solution)
    return [calc, solution]


def read_output_settings(config):
    output = {}
    output['directory'] = config.output.directory
    output['save_data'] = int(config.output.save_data)
    output['save_figures'] = int(config.output.save_figures)
    return output

    
def read_config(filepath):
    # Display status
    sys.stdout.write('\nReading out the configuration file... ')
    # Initialize dictionaries for input parameters 
    exp = {}
    spinA = {}
    spinB = {}
    mode = {}
    simPar = {}
    fitPar = {}
    calc = {}
    solution = {}
    output = {}
    # Read out the config file
    with io.open(filepath) as file:
        config = libconf.load(file)
        # Experimental parameters
        exp = read_experimental_parameters(config)
        # Spin system parameters
        spinA, spinB = read_spin_parameters(config, exp)
        # Calculation mode
        mode = read_calculation_mode(config)
        # Read out the simulation parameters
        if mode['simulation']:
            simPar = read_simulation_parameters(config, exp)
        # Read out fitting parameters 
        if mode['fitting'] or mode['validation']:
            fitPar = read_fitting_parameters(config)
        # Read out calculation settings
        calc, solution = read_calculation_parameters(config, mode, exp, simPar)
        # Output settings
        output = read_output_settings(config)
        # Display status
        sys.stdout.write('[DONE]\n')
    return [exp, spinA, spinB, mode, simPar, fitPar, calc, solution, output]