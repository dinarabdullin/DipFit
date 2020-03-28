'''
calculate_parameter_error.py 

Calculates the confidence interval(s) of the DipFit fitting parameter(s) x and y
	
Optional arguments:
	--threshold         set the threshold in sigma units
    --xopt              set the optimal value of parameter 1
    --yopt              set the optimal value of parameter 2
    --numericalerror    set the numerical error

Requirements: Python3, argparse, wx, numpy, scipy, matplotlib 
'''

import os
import sys
import argparse
import wx
import numpy as np
from scipy.optimize import curve_fit
from functools import partial
import matplotlib.pyplot as plt
import matplotlib.collections as collections
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
rcParams['font.size'] = 24

const = {}
const['variableLabels'] = {
	'r_mean'   : r'$\langle\mathit{r}\rangle$ (nm)',
	'r_width'  : r'$\mathit{\Delta r}$ (nm)', 
	'xi_mean'  : r'$\langle\mathit{\xi}\rangle$ $^\circ$', 
	'xi_width' : r'$\mathit{\Delta\xi}$ $^\circ$', 
	'phi_mean' : r'$\langle\mathit{\varphi}\rangle$ $^\circ$', 
	'phi_width': r'$\mathit{\Delta\varphi}$ $^\circ$', 
	'temp'     : r'Temperature (K)'}


def get_path(message):
    app = wx.App(None) 
    dialog = wx.FileDialog(None, message, wildcard='*.*', style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
    if dialog.ShowModal() == wx.ID_OK:
        path = dialog.GetPath()
    return path
    

def read_validation_file(filepath):
    par1_list = []
    par2_list = []
    score_list = []
    file = open(filepath, 'r')
    next(file)
    for line in file:
        str = line.split()
        ncol = len(str)
        par1_list.append(float(str[0]))
        if ncol == 2:
            score_list.append(float(str[1]))
        elif ncol == 3:
            par2_list.append(float(str[1]))
            score_list.append(float(str[2]))
        else:
            print("Illegible data format!")
    par1_data = np.array(par1_list)
    par2_data = np.array(par2_list)
    score_data = np.array(score_list)  
    return [par1_data, par2_data, score_data]


def gauss_fit(x, a, dx, x0, y0):
    return y0 + a * (1 - np.exp(-0.5 * ((x-x0) / dx)**2))


def calculate_confidence_interval(x_data, y_data, x_name, threshold, x_opt, num_error, filename):
    # Determine the minimal and maximal values of x
    x_min = np.amin(x_data)
    x_max = np.amax(x_data)
    # Create the new x-axis
    Nx = 101
    x_data_interpolated = np.linspace(x_min, x_max, Nx)
    x_data_inc = x_data_interpolated[1] - x_data_interpolated[0]
    x_data_interpolated_lower = x_data_interpolated - x_data_inc * np.ones(x_data_interpolated.size)
    x_data_interpolated_lower[0] = x_data_interpolated[0]
    x_data_interpolated_upper = x_data_interpolated + x_data_inc * np.ones(x_data_interpolated.size)
    x_data_interpolated_upper[-1] = x_data_interpolated[-1]
    # Select the minimal value of y for each x value
    y_data_interpolated = np.zeros(x_data_interpolated.size)
    for i in range(Nx):
        y_data_interpolated[i] = np.amin(y_data[(x_data_interpolated_lower[i] < x_data) & (x_data < x_data_interpolated_upper[i])])  
    # Set the optimal values of x and y
    if np.isnan(x_opt):
        y_opt = np.amin(y_data_interpolated)
        idx_y_opt = np.argmin(y_data_interpolated)
        x_opt = x_data_interpolated[idx_y_opt]
    else:
        idx_x_opt = min(range(len(x_data_interpolated)), key=lambda i: abs(x_data_interpolated[i]-x_opt))
        y_opt = y_data_interpolated[idx_x_opt]
    # Fit the dependence y(x) using the parameters a and x0 of the function gauss_fit:
    func = partial(gauss_fit, x0=x_opt, y0=y_opt)
    popt, pcov = curve_fit(func, x_data_interpolated, y_data_interpolated)
    # Determine the confidence interval
    a = abs(popt[0])
    dx_opt = abs(popt[1])
    x_error = threshold * dx_opt
    x_left = x_opt - x_error
    x_right = x_opt + x_error
    if a <= 0.2 * np.amin(y_data_interpolated):
        x_left = x_min
        x_right = x_max
        x_error_final = np.nan
        y_dev = np.nan
    else:
        if (x_left >= x_min) and (x_right <= x_max):
            x_error_final = x_error
            idx_left = min(range(len(x_data_interpolated)), key=lambda i: abs(x_data_interpolated[i]-x_left))
            idx_right = min(range(len(x_data_interpolated)), key=lambda i: abs(x_data_interpolated[i]-x_right))
            y_left = y_data_interpolated[idx_left]
            y_right = y_data_interpolated[idx_right]
            y_dev_left = abs(y_opt - y_left)
            y_dev_right = abs(y_opt - y_right)
            y_dev = np.amin([y_dev_left, y_dev_right])
        elif (x_left < x_min) and (x_right <= x_max):
            x_left = x_min
            x_error_final = x_error
            idx_right = min(range(len(x_data_interpolated)), key=lambda i: abs(x_data_interpolated[i]-x_right))
            y_right = y_data_interpolated[idx_right]
            y_dev = abs(y_opt - y_right)
        elif (x_left >= x_min) and (x_right > x_max):
            x_right = x_max
            x_error_final = x_error
            idx_left = min(range(len(x_data_interpolated)), key=lambda i: abs(x_data_interpolated[i]-x_left))
            y_left = y_data_interpolated[idx_left]
            y_dev = abs(y_opt - y_left)
        elif (x_left < x_min) and (x_right > x_max):
            x_left = x_min
            x_right = x_max
            x_error_final = np.nan
            y_dev = np.nan
    print('Confidence interval of %s: %f (%d*sigma level)' % (x_name,x_error_final, threshold))
    # Check that the numerical error is inside the RMSD threshold
    if not np.isnan(num_error) and not np.isnan(y_dev):
        if num_error > y_dev:
            sys.stdout.write('Warning: The RMSD threshold is under the numerical error! Increase the threshold.\n')
    # Plot the results  
    x_fit = x_data_interpolated
    y_fit = func(x_data_interpolated, *popt)
    plot_confidence_interval(x_data, y_data, x_data_interpolated, y_data_interpolated, x_fit, y_fit, [x_left, x_right], x_name, filename)
    return x_error_final


def plot_confidence_interval(x, y, xi, yi, xf, yf, x_ci, x_name, filename):    
    # Set the maximum and the minimum of y
    ymin = np.nanmin(y)
    ymax = 2 * ymin
    # Plot the figure
    fig = plt.figure()
    fig.patch.set_facecolor('white')
    ax = plt.gca()
    ax.scatter(x, y, c=y, cmap='jet_r', vmin=ymin, vmax=ymax, label='data')
    #ax.plot(x, y, 'ro', label='data')
    #ax.plot(xi, yi, 'b-', label='interpolated data')
    ax.plot(xf, yf, 'k-', label='Gaussian fit')
    ax.axvspan(x_ci[0], x_ci[1], facecolor='gray', alpha=0.4, label='confidence interval')
    ax.set_xlim(round(np.amin(x),1), round(np.amax(x),1))
    plt.xlabel(const['variableLabels'][x_name])
    plt.ylabel('RMSD')
    #plt.legend()
    plt.tight_layout()
    plt.draw()
    plt.show(block=False)
    plt.savefig(filename, format='png', dpi=600)


def keep_figures_live():
	plt.show()	


if __name__ == '__main__':
    # Input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-th', '--threshold', type=float, help="set the threshhold in sigma units")
    parser.add_argument('-xo', '--xopt', type=float, help="set the optimal value of parameter 1")
    parser.add_argument('-yo', '--yopt', type=float, help="set the optimal value of parameter 2")
    parser.add_argument('-ne', '--numericalerror', type=float, help="set the numerical error")
    args = parser.parse_args()
    threshold = 1
    if args.threshold:
        threshold = args.threshold
    par1_opt = np.nan
    if args.xopt:
        par1_opt = args.xopt
    par2_opt = np.nan
    if args.yopt:
        par2_opt = args.yopt
    num_error = np.nan
    if args.numericalerror:
        num_error = args.numericalerror   
    # Read data
    filepath = get_path("Open the validation file...")
    filedir = os.path.dirname(filepath)
    filename = os.path.basename(filepath[:-4])
    filename_parts = filename.split("-")
    par1, par2, score = read_validation_file(filepath)
    # Calculate the confidence interval of parameter 1
    par1_name = filename_parts[1]
    filename = filedir + "/confidence_interval_" + par1_name + ".png"
    par1_error = calculate_confidence_interval(par1, score, par1_name, threshold, par1_opt, num_error, filename)
    # Calculate the confidence interval of parameter 2
    if not (par2.size == 0):
        par2_name = filename_parts[2]
        filename = filedir + "/confidence_interval_" + par2_name + ".png"
        par2_error = calculate_confidence_interval(par2, score, par2_name, threshold, par2_opt, num_error, filename)
    print('Done!')
    keep_figures_live()    