'''
plot_error_analysis_data.py 

Plots the validation data of DipFit
	
Optional arguments:
	--fontsize          set font size
    --threshold         set the threshold
    --numerical_error   set the numerical error
    --rmsd_min          set the minimal RMSD

Requirements: Python3, argparse, wx, numpy, scipy, matplotlib 
'''

import argparse
import os
import io
import sys
import wx
import numpy as np
from copy import deepcopy
from scipy.interpolate import griddata
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['xtick.direction'] = 'out'
rcParams['ytick.direction'] = 'out'
rcParams['axes.facecolor']= 'white'
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = 'Arial'
rcParams['lines.linewidth'] = 2
rcParams['xtick.major.size'] = 8
rcParams['xtick.major.width'] = 1.5
rcParams['ytick.major.size'] = 8
rcParams['ytick.major.width'] = 1.5
rcParams['font.size'] = 18


error_analysis_files = [
    {'filename': 'parameter_errors-r_mean-r_width.dat',     '2d': True, 'x1': 'r_mean',     'x2': 'r_width'     },
    {'filename': 'parameter_errors-xi_mean-xi_width.dat',   '2d': True, 'x1': 'xi_mean',    'x2': 'xi_width'    },
    {'filename': 'parameter_errors-phi_mean-phi_width.dat', '2d': True, 'x1': 'phi_mean',   'x2': 'phi_width'   },
    {'filename': 'parameter_errors-xi_mean-phi_mean.dat',   '2d': True, 'x1': 'xi_mean',    'x2': 'phi_mean'    },
    {'filename': 'parameter_errors-r_width-xi_mean.dat',    '2d': True, 'x1': 'r_width',    'x2': 'xi_mean'     },
    {'filename': 'parameter_errors-r_width-xi_width.dat',   '2d': True, 'x1': 'r_width',    'x2': 'xi_width'    },
    {'filename': 'parameter_errors-r_width-phi_mean.dat',   '2d': True, 'x1': 'r_width',    'x2': 'phi_mean'    },
    {'filename': 'parameter_errors-r_width-phi_width.dat',  '2d': True, 'x1': 'r_width',    'x2': 'phi_width'   },
    {'filename': 'parameter_errors-r_mean.dat',             '2d': False,'x1': 'r_mean',     'x2': ''            },
    {'filename': 'parameter_errors-r_width.dat',            '2d': False,'x1': 'r_width',    'x2': ''            },
    {'filename': 'parameter_errors-xi_mean.dat',            '2d': False,'x1': 'xi_mean',    'x2': ''            },
    {'filename': 'parameter_errors-xi_width.dat',           '2d': False,'x1': 'xi_width',   'x2': ''            },
    {'filename': 'parameter_errors-phi_mean.dat',           '2d': False,'x1': 'phi_mean',   'x2': ''            },
    {'filename': 'parameter_errors-phi_width.dat',          '2d': False,'x1': 'phi_width',  'x2': ''            },
    {'filename': 'parameter_errors-temp.dat',               '2d': False,'x1': 'temp',       'x2': ''            },
]


const = {}
const['variable_names'] = [
	'r_mean',
	'r_width', 
	'xi_mean', 
	'xi_width', 
	'phi_mean', 
	'phi_width', 
	'temp']
const['long_variable_names'] = {
	'r_mean'   : 'r mean (nm)',
	'r_width'  : 'r width (nm)', 
	'xi_mean'  : 'xi mean (deg)',
	'xi_width' : 'xi width (deg)', 
	'phi_mean' : 'phi mean (deg)', 
	'phi_width': 'phi mean (deg)',
	'temp'     : 'temperature (K)'}
const['variable_scales'] = {
	'r_mean'   : 1.0,
	'r_width'  : 1.0, 
	'xi_mean'  : np.pi / 180.0, 
	'xi_width' : np.pi / 180.0, 
	'phi_mean' : np.pi / 180.0, 
	'phi_width': np.pi / 180.0, 
	'temp'     : 1.0}
const['variable_labels'] = {
	'r_mean'   : r'$\langle\mathit{r}\rangle$ (nm)',
	'r_width'  : r'$\mathit{\Delta r}$ (nm)', 
	'xi_mean'  : r'$\langle\mathit{\xi}\rangle$ $^\circ$', 
	'xi_width' : r'$\mathit{\Delta\xi}$ $^\circ$', 
	'phi_mean' : r'$\langle\mathit{\varphi}\rangle$ $^\circ$', 
	'phi_width': r'$\mathit{\Delta\varphi}$ $^\circ$', 
	'temp'     : r'Temperature (K)'}


# No. of plots: 1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16
alignement = [[1,1], [1,2], [1,3], [2,2], [2,3], [2,3], [2,4], [2,4], [3,3], [3,4], [3,4], [3,4], [4,4], [4,4], [4,4], [4,4]]


def get_path(message):
    app = wx.App(None) 
    dialog = wx.FileDialog(None, message, wildcard='*.*', style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
    if dialog.ShowModal() == wx.ID_OK:
        path = dialog.GetPath()
    return path


def load_error_analysis_data(directory):
    nfiles = len(error_analysis_files)
    parameters = []
    score_vs_parametes = []
    for i in range(nfiles):
        filename = directory + error_analysis_files[i]['filename']
        # check if the file exists
        if os.path.isfile(filename):
            # set parameters
            if error_analysis_files[i]['2d'] == False:
                variables = [error_analysis_files[i]['x1']]
            else:
                variables = [error_analysis_files[i]['x1'], error_analysis_files[i]['x2']]
            parameters.append(variables)
            # read the data from the file
            data = np.genfromtxt(filename, skip_header=1)
            nx = data.shape[0]
            score_vs_varibales = {}
            if error_analysis_files[i]['2d'] == False:
                name1 = error_analysis_files[i]['x1']
                score_vs_varibales[name1] = [data[j][0] * const['variable_scales'][name1] for j in range(nx)]
                score_vs_varibales['score'] = [data[j][1] for j in range(nx)]
            else:
                name1 = error_analysis_files[i]['x1']
                name2 = error_analysis_files[i]['x2']
                score_vs_varibales[name1] = [data[j][0] * const['variable_scales'][name1] for j in range(nx)]
                score_vs_varibales[name2] = [data[j][1] * const['variable_scales'][name2] for j in range(nx)]
                score_vs_varibales['score'] = [data[j][2] for j in range(nx)]
            score_vs_parametes.append(score_vs_varibales)
    return [parameters, score_vs_parametes]


def chunkstring(string, length):
    return (string[0+i:length+i] for i in range(0, len(string), length))


def load_optimized_parameters(filepath):
    optimized_parameters = {}
    count = 0
    file = open(filepath, 'r')
    # Skip a header
    next(file)
    for line in file:
        str = list(chunkstring(line, 16))
        parameter = {}
        name = const['variable_names'][count] 
        parameter['longname'] = str[0].strip()
        parameter['value'] = float(str[1]) * const['variable_scales'][name]
        parameter['optimized'] = str[2].strip()
        parameter['precision'] = float(str[3]) * const['variable_scales'][name]
        optimized_parameters[name] = parameter
        count += 1
    return optimized_parameters


def calculate_best_score(score_min, parameters, score_vs_parameters):
    best_score = 0
    if score_min:
        best_score = score_min
    else:
        best_scores = []
        for i in range(len(parameters)):
            best_scores.append(np.amin(score_vs_parameters[i]['score']))
        best_scores = np.array(best_scores) 
        best_score = np.amin(best_scores)
    return best_score


def calculate_score_threshold(best_score, threshold, numerical_error):
    score_threshold = best_score * threshold
    if (numerical_error > (score_threshold - best_score)):
        score_threshold = best_score + numerical_error
        sys.stdout.write('Warning: The RMSD threshold is below the numerical error! The RMSD threshold will be increased to the numerical error level.\n')
    sys.stdout.write("RMSD threshold = %f\n" % (score_threshold))
    return score_threshold


def calculate_parameter_errors(variables, score_vs_parameters, best_parameters, score_threshold):
    parameter_errors = {}
    Ne = len(score_vs_parameters)
    for i in range(Ne):
        for name in variables[i]:
            variable_values = np.array(score_vs_parameters[i][name])
            score_values = np.array(score_vs_parameters[i]['score'])
            best_parameter = best_parameters[name]['value']
            parameter_error = calculate_parameter_error(variable_values, score_values, best_parameter, score_threshold)
            if name in parameter_errors:
                if not np.isnan(parameter_error) and not np.isnan(parameter_errors[name]):
                    if (parameter_error > parameter_errors[name]):
                        parameter_errors[name] = parameter_error
            else:
                parameter_errors[name] = parameter_error
    return parameter_errors


def calculate_parameter_error(x_data, y_data, x_opt, threshold):
    Ns = x_data.size
    # Determine the minimal and maximal values of x
    x_min = np.amin(x_data)
    x_max = np.amax(x_data)
    # Set the optimal values of x and y
    if np.isnan(x_opt):
        y_opt = np.amin(y_data)
        idx_y_opt = np.argmin(y_data)
        x_opt = x_data[idx_y_opt]
    else:
        idx_x_opt = min(range(len(x_data)), key=lambda i: abs(x_data[i]-x_opt))
        y_opt = y_data[idx_x_opt]
    # Sort x in ascending order
    x_sorted, y_sorted = zip(*sorted(zip(x_data, y_data)))
    # Determine the uncertainty ranges of x
    idx_x_selected = []
    for i in range(Ns):
        if y_sorted[i] < threshold:
            idx_x_selected.append(i)
    x_left = x_sorted[idx_x_selected[0]]
    x_right = x_sorted[idx_x_selected[-1]]
    # Determine the error of x_opt
    x_dev_left = abs(x_opt - x_left)
    x_dev_right = abs(x_opt - x_right)
    x_error = np.amax([x_dev_left, x_dev_right])
    if (x_left == x_min) and (x_right == x_max):
        x_error = np.nan
    return x_error


def include_errors(optimized_parameters, parameter_errors):
    optimized_parameters_with_errors = deepcopy(optimized_parameters)
    for name in parameter_errors:
        if not np.isnan(parameter_errors[name]):
            optimized_parameters_with_errors[name]['precision'] = parameter_errors[name]
    return optimized_parameters_with_errors


def print_optimized_parameters(optimized_parameters):
        sys.stdout.write('Optimized fitting parameters:\n')
        sys.stdout.write("{0:<16s} {1:<16s} {2:<16s} {3:<16s}\n".format('Parameter', 'Value', 'Optimized', 'Precision (+/-)'))
        for name in const['variable_names']:
            parameter = optimized_parameters[name]
            sys.stdout.write("{0:<16s} ".format(parameter['longname']))
            sys.stdout.write("{0:<16.3f} ".format(parameter['value'] / const['variable_scales'][name]))
            sys.stdout.write("{0:<16s} ".format(parameter['optimized']))
            sys.stdout.write("{0:<16.3f} \n".format(parameter['precision'] / const['variable_scales'][name]))
        sys.stdout.write('\n')


def plot_score_vs_parameters(variables, score_vs_parameters, score_threshold, best_parameters, save_figure=False, filename=''): 
    Ne = len(variables)
    c = 1
    fig = plt.figure(figsize=(18,9), facecolor='w', edgecolor='w')
    for i in range(Ne):
        dim = len(variables[i])
        if (dim == 1):
            plt.subplot(alignement[Ne-1][0], alignement[Ne-1][1], c)
            c = c + 1
            im = plot_1d(fig, variables[i], score_vs_parameters[i], score_threshold, best_parameters)
        elif (dim == 2):
            plt.subplot(alignement[Ne-1][0], alignement[Ne-1][1], c)
            c = c + 1
            im = plot_2d(fig, variables[i], score_vs_parameters[i], score_threshold, best_parameters)
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15, right=0.80, top=0.9)
    cax = plt.axes([0.85, 0.3, 0.02, 0.4]) # left, bottom, width, height  
    plt.colorbar(im, cax=cax, orientation='vertical')
    plt.text(0.90, 1.05, 'RMSD', transform=cax.transAxes)
    plt.draw()
    plt.show(block=False)
    if save_figure:
        plt.savefig(filename, format='png', dpi=600)


def plot_1d(fig, variables, score_vs_parameters, score_threshold, best_parameters):
    # Read out the values of fitting parameters and RMSD values
    name1 = variables[0]
    x1 = [(x / const['variable_scales'][name1]) for x in score_vs_parameters[name1]]
    y = score_vs_parameters['score']
    # Read out the optimized values of fitting parameters
    x1_opt = best_parameters[name1]['value'] / const['variable_scales'][name1]
    y_opt = np.amin(y)
    # Set the maximum and the minimum of y
    ymin = score_threshold
    ymax = 2 * score_threshold
    # Plot the figure
    axes = fig.gca()
    im = axes.scatter(x1, y, c=y, cmap='jet_r', vmin=ymin, vmax=ymax)
    axes.set_xlim(round(np.amin(x1),1), round(np.amax(x1),1))
    axes.set_xlabel(const['variable_labels'][name1])
    axes.set_ylabel('RMSD')
    axes.plot(x1_opt, y_opt, color='black', marker='o', markerfacecolor='white', markersize=12, clip_on=False)
    plt.margins(0.0)
    # Make the axes of equal scale
    xl, xh = axes.get_xlim()
    yl, yh = axes.get_ylim()
    axes.set_aspect( (xh-xl)/(yh-yl) )
    return im


def plot_2d(fig, variables, score_vs_parameters, score_threshold, best_parameters):
    # Read out the values of fitting parameters and RMSD values
    name1 = variables[0]
    name2 = variables[1]
    x1 = [(x / const['variable_scales'][name1]) for x in score_vs_parameters[name1]]	
    x2 = [(x / const['variable_scales'][name2]) for x in score_vs_parameters[name2]]	
    y = score_vs_parameters['score']
    # Read out th optimized values of fitting parameters
    x1_opt = best_parameters[name1]['value'] / const['variable_scales'][name1]
    x2_opt = best_parameters[name2]['value'] / const['variable_scales'][name2]
    # Interpolate the data on a regular grid
    x1min = np.min(x1)
    x1max = np.max(x1)
    x2min = np.min(x2)
    x2max = np.max(x2)
    x1r = np.linspace(x1min, x1max, num=200) 
    x2r = np.linspace(x2min, x2max, num=200)
    X, Y = np.mgrid[x1min:x1max:200j, x2min:x2max:200j]
    Z = griddata((x1, x2), y, (X, Y), method='linear')
    # Set the maximum and the minimum of Z
    zmin = score_threshold
    zmax = 2 * score_threshold
    # Plot the figure
    axes = fig.gca()
    im = axes.pcolor(X, Y, Z, cmap='jet_r', vmin=zmin, vmax=zmax)
    axes.set_xlim(np.amin(x1), np.amax(x1))
    axes.set_ylim(np.amin(x2), np.amax(x2))
    axes.set_xlabel(const['variable_labels'][name1])
    axes.set_ylabel(const['variable_labels'][name2])
    axes.plot(x1_opt, x2_opt, color='black', marker='o', markerfacecolor='white', markersize=12, clip_on=False)
    plt.margins(0.0)
    # Make the axis of equal scale
    xl, xh = axes.get_xlim()
    yl, yh = axes.get_ylim()
    axes.set_aspect( (xh-xl)/(yh-yl) )
    # Make ticks
    plt.xticks(np.linspace(round(x1min,1), round(x1max,1), 3))
    plt.yticks(np.linspace(round(x2min,1), round(x2max,1), 3))
    axes.xaxis.labelpad = 0
    axes.yaxis.labelpad = 0
    return im


def keep_figures_live():
	plt.show()	

	
if __name__ == '__main__':
    # Input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-fs', '--fontsize', type=int, help="set the font size")
    parser.add_argument('-th', '--threshold', type=float, help="set the threshhold")
    parser.add_argument('-ne', '--numerical_error', type=float, help="set the numerical error")
    parser.add_argument('-rm', '--rmsd_min', type=float, help="set the minimal RMSD")
    args = parser.parse_args()
    if args.fontsize:
        fontsize = args.fontsize
        matplotlib.rcParams.update({'font.size': fontsize}) 
    threshold = 1.2
    if args.threshold:
        threshold = 0.01 * args.threshold
    numerical_error = 0
    if args.numerical_error:
        numerical_error = args.numerical_error 
    score_min = 0
    if args.rmsd_min:
        score_min = args.rmsd_min    
    # Read the results of the error analysis
    filepath = get_path("Open file with the the results of error analysis...")
    directory = os.path.dirname(filepath) + '/'
    variables, score_vs_parameters = load_error_analysis_data(directory)
    # Read the optimized fitting parameters
    filepath2 = get_path("Open file with the optimized fitting parameters...")
    best_parameters = load_optimized_parameters(filepath2)
    # Calculate the best score
    best_score = calculate_best_score(score_min, variables, score_vs_parameters)
    # Calculate the threshold
    score_threshold = calculate_score_threshold(best_score, threshold, numerical_error)
    # Calculate the errors of the fitting parameters
    parameter_errors = calculate_parameter_errors(variables, score_vs_parameters, best_parameters, score_threshold)
    # Update the optimized fitting parameters
    best_parameters_with_errors = include_errors(best_parameters, parameter_errors)
    # Display the optmized fitting parameters
    print_optimized_parameters(best_parameters_with_errors)
    # Plot the validation data
    filename = directory + 'parameter_errors_formated.png'
    plot_score_vs_parameters(variables, score_vs_parameters, score_threshold, best_parameters, True, filename)	
    keep_figures_live()