'''
Read the input data from the config file
'''

import io
import sys
import libconf
import numpy as np
from input.load_spectrum import load_spectrum
from input.load_timetrace import load_timetrace
from input.load_optimized_parameters import load_optimized_parameters
from input.symmetric_boundaries import symmetric_boundaries
from mathematics.find_max_on_interval import find_max_on_interval
from supplement.constants import const


def read_calculation_mode(config):
    mode = {}
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
    else:
        raise ValueError('Illelgible value of mode!')
        sys.exit(1)
    return mode


def read_experimental_parameters(config):
    expData = {}  
    expData['path_spectrum'] = config.path_spectrum
    expData['path_timetrace'] = config.path_timetrace
    expData['t'] = []
    expData['sig'] = []
    expData['f'] = []
    expData['spc'] = []
    if not (expData['path_spectrum']==""):
        f, spc = load_spectrum(expData['path_spectrum'])
        expData['f'], expData['spc'] = symmetric_boundaries(f, spc)
    if not (expData['path_timetrace']==""):
        expData['t'], expData['sig'] = load_timetrace(expData['path_timetrace'])
    return expData


def read_spin_parameters(config):
    spinA = {}
    spinB = {}
    spinA['type'] = config.spinA.type
    spinA['g'] = np.array([float(config.spinA.g[0]), 
                           float(config.spinA.g[1]), 
                           float(config.spinA.g[2])])
    spinB['type'] = config.spinB.type
    spinB['g'] = np.array([float(config.spinB.g[0]), 
                           float(config.spinB.g[1]), 
                           float(config.spinB.g[2])])
    if (spinA['type'] == "anisotropic") and (spinB['type'] == "isotropic"):
        raise ValueError('Only SpinB can be anisotropic!')
        sys.exit(1)
    if (spinA['type'] == "anisotropic") and (spinB['type'] == "anisotropic"):
        raise ValueError('The case of two anisotropic spins is not supported yet!')
        sys.exit(1)
    return [spinA, spinB]


def read_simulation_settings(config, expData):
    simSettings = {}
    simSettings['modes'] = {}
    simSettings['modes']['spc'] = int(config.simulation_modes.spc)
    simSettings['modes']['timetrace'] = int(config.simulation_modes.timetrace)
    simSettings['modes']['spc_vs_theta'] = int(config.simulation_modes.spc_vs_theta) 
    simSettings['modes']['spc_vs_xi'] = int(config.simulation_modes.spc_vs_xi)
    simSettings['modes']['spc_vs_phi'] = int(config.simulation_modes.spc_vs_phi)
    simSettings['modes']['spc_vs_temp'] = int(config.simulation_modes.spc_vs_temp)
    simSettings['variables'] = {}
    simSettings['variables']['r_mean'] = float(config.simulation_variables.r_mean)
    simSettings['variables']['r_width'] = float(config.simulation_variables.r_width)
    simSettings['variables']['xi_mean'] = const['deg2rad'] * float(config.simulation_variables.xi_mean)
    simSettings['variables']['xi_width'] = const['deg2rad'] * float(config.simulation_variables.xi_width)
    simSettings['variables']['phi_mean'] = const['deg2rad'] * float(config.simulation_variables.phi_mean)
    simSettings['variables']['phi_width'] = const['deg2rad'] * float(config.simulation_variables.phi_width)
    simSettings['variables']['temp'] = float(config.simulation_variables.temp)
    simSettings['settings'] = {}
    simSettings['settings']['theta_bins'] = []
    simSettings['settings']['xi_bins'] = []
    simSettings['settings']['phi_bins'] = []
    simSettings['settings']['temp_bins'] = []
    if simSettings['modes']['spc_vs_theta']:
        simSettings['settings']['theta_bins'] = np.linspace(float(config.simulation_settings.theta_ranges[0]), 
                                                            float(config.simulation_settings.theta_ranges[1]), 
                                                            int(config.simulation_settings.theta_ranges[2]))
    if simSettings['modes']['spc_vs_xi']:
        simSettings['settings']['xi_bins'] = np.linspace(float(config.simulation_settings.xi_ranges[0]), 
                                                         float(config.simulation_settings.xi_ranges[1]), 
                                                         int(config.simulation_settings.xi_ranges[2]))
    if simSettings['modes']['spc_vs_phi']:
        simSettings['settings']['phi_bins'] = np.linspace(float(config.simulation_settings.phi_ranges[0]), 
                                                          float(config.simulation_settings.phi_ranges[1]), 
                                                          int(config.simulation_settings.phi_ranges[2]))
    if simSettings['modes']['spc_vs_temp']:
        simSettings['settings']['temp_bins'] = np.linspace(float(config.simulation_settings.temp_ranges[0]), 
                                                           float(config.simulation_settings.temp_ranges[1]), 
                                                           int(config.simulation_settings.temp_ranges[2]))
    if (expData['f'] == []):
        simSettings['settings']['faxis_normalized'] = int(config.simulation_settings.faxis_normalized)
    else:
        simSettings['settings']['faxis_normalized'] = 0
    simSettings['settings']['plot_3d'] = int(config.simulation_settings.plot_3d)
    return simSettings


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


def read_fitting_settings(config, expData): 
    fitSettings = {}
    fitSettings['variables'] = {}
    fitSettings['variables']['indices'] = {}
    fitSettings['variables']['bounds'] = []
    fitSettings['variables']['fixed'] = {}
    fitSettings['variables']['indices'], fitSettings['variables']['bounds'], fitSettings['variables']['fixed'], fitSettings['variables']['size'] = read_fitting_variables(config.fitting_variables)
    fitSettings['settings'] = {}
    fitSettings['settings']['fitted_data'] = config.fitting_settings.fitted_data
    fitSettings['settings']['display_graphics'] = int(config.fitting_settings.display_graphics)
    fitSettings['settings']['method'] = config.fitting_settings.method
    if (fitSettings['settings']['method'] == "genetic"):
        fitSettings['settings']['num_generations'] = int(config.fitting_settings.num_generations)
        fitSettings['settings']['size_generation'] = int(config.fitting_settings.size_generation)
        fitSettings['settings']['prob_crossover'] = float(config.fitting_settings.prob_crossover)
        fitSettings['settings']['prob_mutation'] = float(config.fitting_settings.prob_mutation)
    # Check the spin system
    if (fitSettings['settings']['fitted_data'] == 'spectrum') and (expData['f'] == []):
        raise ValueError('No spectrum was specified!')
        sys.exit(1)
    elif (fitSettings['settings']['fitted_data'] == 'timetrace') and (expData['t'] == []):
        raise ValueError('No time trace was specified!')
        sys.exit(1)
    return fitSettings


def read_validation_settings(config, mode):
    valSettings = {}
    valSettings['variables'] = config.validation.variables
    valSettings['Ns'] = int(config.validation.Ns)
    valSettings['threshold'] = float(config.validation.threshold)
    path_optimized_parameters = config.validation.path_optimized_parameters
    if not (path_optimized_parameters == "") and (mode['validation']):
        valSettings['optimized_parameters'] = load_optimized_parameters(path_optimized_parameters)
    return valSettings


def read_calculation_settings(config, expData):
    calcSettings = {}
    calcSettings['Ns'] = int(config.calculation_settings.Ns) 
    calcSettings['r_distr'] = config.calculation_settings.r_distr
    calcSettings['xi_distr'] = config.calculation_settings.xi_distr
    calcSettings['phi_distr'] = config.calculation_settings.phi_distr
    calcSettings['f_min'] = float(config.calculation_settings.fmin)
    calcSettings['f_max'] = float(config.calculation_settings.fmax)
    calcSettings['t_min'] = float(config.calculation_settings.tmin)
    calcSettings['t_max'] = float(config.calculation_settings.tmax)
    calcSettings['g_selectivity'] = int(config.calculation_settings.g_selectivity)
    if calcSettings['g_selectivity']:
        calcSettings['magnetic_field'] = float(config.calculation_settings.magnetic_field)
    # Set the maximal amplitude of the spectrum and the maximal dipolar frequency
    if not (expData['f'] == []):
        if calcSettings['f_max']:
            calcSettings['spc_max'] = find_max_on_interval(expData['spc'], expData['f'], calcSettings['f_min'], calcSettings['f_max'])   
        else:
            calcSettings['f_max'] = np.amax(expData['f'])
            calcSettings['spc_max'] = find_max_on_interval(expData['spc'], expData['f'], calcSettings['f_min'], calcSettings['f_max'])       
    else:
        calcSettings['spc_max'] = 1.0
    # Set the maximal time point    
    if not (expData['t'] == []):
        if not calcSettings['t_max']:
            calcSettings['t_max'] = np.amax(expData['t'])
    return calcSettings


def read_output_settings(config):
    outputSettings = {}
    outputSettings['directory'] = config.output.directory
    outputSettings['save_data'] = int(config.output.save_data)
    outputSettings['save_figures'] = int(config.output.save_figures)
    return outputSettings

    
def read_config(filepath):
    sys.stdout.write('\n') 
    sys.stdout.write('Reading out the configuration file... ') 
    mode = {}
    expData = {}
    spinA = {}
    spinB = {}
    simSettings = {}
    fitSettings = {}
    valSettings = {}
    calcSettings = {}
    outputSettings = {}
    with io.open(filepath) as file:
        config = libconf.load(file)
        mode = read_calculation_mode(config)
        expData = read_experimental_parameters(config)
        spinA, spinB = read_spin_parameters(config)
        if mode['simulation']:
            simSettings = read_simulation_settings(config, expData)
        if mode['fitting'] or mode['validation']:
            fitSettings = read_fitting_settings(config, expData)
            valSettings = read_validation_settings(config, mode)
        calcSettings = read_calculation_settings(config, expData)
        outputSettings = read_output_settings(config)
    sys.stdout.write('[DONE]\n\n')
    return [mode, expData, spinA, spinB, simSettings, fitSettings, valSettings, calcSettings, outputSettings]
    