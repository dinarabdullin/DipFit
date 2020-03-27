'''
plot_validation_data.py 

Plots the validation data of DipFit
	
Optional arguments:
	--fontsize      set font size

Requirements: Python3, argparse, wx, numpy, scipy, matplotlib 
'''

import argparse
import os
import io
import sys
import wx
import numpy as np
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


validation_files = [
    {'filename': 'parameter_errors-r_mean-r_width.dat',     '2d': True, 'x1': 'r_mean',     'x2': 'r_width'     },
    {'filename': 'parameter_errors-xi_mean-xi_width.dat',   '2d': True, 'x1': 'xi_mean',    'x2': 'xi_width'    },
    {'filename': 'parameter_errors-phi_mean-phi_width.dat', '2d': True, 'x1': 'phi_mean',   'x2': 'phi_width'   },
    {'filename': 'parameter_errors-xi_mean-phi_mean.dat',   '2d': True, 'x1': 'xi_mean',    'x2': 'phi_mean'    },
    {'filename': 'parameter_errors-xi_mean-r_width.dat',    '2d': True, 'x1': 'xi_mean',    'x2': 'r_width'     },
    {'filename': 'parameter_errors-xi_width-r_width.dat',   '2d': True, 'x1': 'xi_width',   'x2': 'r_width'     },
    {'filename': 'parameter_errors-phi_mean-r_width.dat',   '2d': True, 'x1': 'phi_mean',   'x2': 'r_width'     },
    {'filename': 'parameter_errors-phi_width-r_width.dat',  '2d': True, 'x1': 'phi_width',  'x2': 'r_width'     },
    {'filename': 'parameter_errors-temp.dat',               '2d': False,'x1': 'temp',       'x2': ''            },
]


# No. of plots: 1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16
alignement = [[1,1], [1,2], [1,3], [2,2], [2,3], [2,3], [2,4], [2,4], [3,3], [3,4], [3,4], [3,4], [4,4], [4,4], [4,4], [4,4]]


const= {}
const['variable_labels'] = {
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


def read_validation_data(directory):
    nfiles = len(validation_files)
    par = []
    score_vs_par = []
    for i in range(nfiles):
        filename = directory + validation_files[i]['filename']
        # check if the file exists
        if os.path.isfile(filename):
            # set parameters
            if validation_files[i]['2d'] == False:
                single_par = [validation_files[i]['x1']]
            else:
                single_par = [validation_files[i]['x1'], validation_files[i]['x2']]
            par.append(single_par)
            # read the data from the file
            data = np.genfromtxt(filename, skip_header=1)
            nx = data.shape[0]
            single_score_vs_par = {}
            if validation_files[i]['2d'] == False:
                name1 = validation_files[i]['x1']
                single_score_vs_par[name1] = [data[j][0] for j in range(nx)]
                single_score_vs_par['score'] = [data[j][1] for j in range(nx)]
            else:
                name1 = validation_files[i]['x1']
                name2 = validation_files[i]['x2']
                single_score_vs_par[name1] = [data[j][0] for j in range(nx)]
                single_score_vs_par[name2] = [data[j][1] for j in range(nx)]
                single_score_vs_par['score'] = [data[j][2] for j in range(nx)]
            score_vs_par.append(single_score_vs_par)
    return [par, score_vs_par]


def calculate_best_score(par, score_vs_par):
    best_scores = []
    for i in range(len(par)):
        best_scores.append(np.amin(score_vs_par[i]['score']))
    best_scores = np.array(best_scores) 
    best_score = np.amin(best_scores)
    return best_score


def plot_1d(fig, score_vs_par, var):
    # Read out the values of a gene and the corresponding score values
    name1 = var[0]
    x1 = score_vs_par[name1]
    y = score_vs_par['score']
    # Set the maximum and the minimum of y
    ymin = np.nanmin(y)
    ymax = 2 * ymin
    # Plot the figure
    axes = fig.gca()
    im = axes.scatter(x1, y, c=y, cmap='jet_r', vmin=ymin, vmax=ymax)
    axes.set_xlim(round(np.amin(x1),1), round(np.amax(x1),1))
    axes.set_xlabel(const['variable_labels'][name1])
    axes.set_ylabel('RMSD')
    # Make the axes of equal scale
    xl, xh = axes.get_xlim()
    yl, yh = axes.get_ylim()
    axes.set_aspect( (xh-xl)/(yh-yl) )
    return im


def plot_2d(fig, score_vs_par, var):
    # Read out the values of variables
    name1 = var[0]
    name2 = var[1]
    x1 = score_vs_par[name1]	
    x2 = score_vs_par[name2]	
    y = score_vs_par['score']
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
    zmin = np.nanmin(Z)
    zmax = 2 * zmin
    # Plot the figure
    axes = fig.gca()
    im = axes.pcolor(X, Y, Z, cmap='jet_r', vmin=zmin, vmax=zmax)
    axes.set_xlim(np.amin(x1), np.amax(x1))
    axes.set_ylim(np.amin(x2), np.amax(x2))
    axes.set_xlabel(const['variable_labels'][name1])
    axes.set_ylabel(const['variable_labels'][name2])
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


def plot_score_vs_par(score_vs_par, par, save_figure=False, filename=''): 
    M = len(par)
    N = 1
    fig = plt.figure(figsize=(18,9), facecolor='w', edgecolor='w')
    for i in range(M):
        dim = len(par[i])
        if (dim == 1):
            plt.subplot(alignement[M-1][0], alignement[M-1][1], N)
            N = N + 1
            im = plot_1d(fig, score_vs_par[i], par[i])
        elif (dim == 2):
            plt.subplot(alignement[M-1][0], alignement[M-1][1], N)
            N = N + 1
            im = plot_2d(fig, score_vs_par[i], par[i])
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15, right=0.80, top=0.9)
    cax = plt.axes([0.85, 0.3, 0.02, 0.4]) # left, bottom, width, height  
    plt.colorbar(im, cax=cax, orientation='vertical')
    plt.text(0.90, 1.05, 'RMSD', transform=cax.transAxes)
    plt.draw()
    plt.show(block=False)
    if save_figure:
        plt.savefig(filename, format='png', dpi=600)


def keep_figures_live():
	plt.show()	

	
if __name__ == '__main__':
    # Input arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-fs', '--fontsize', type=int, help="set the font size")
    args = parser.parse_args()
    if args.fontsize:
        fontsize = args.fontsize
        matplotlib.rcParams.update({'font.size': fontsize}) 
    # Read the validation files
    filepath = get_path("Open the directory with validation data...")
    directory = os.path.dirname(filepath) + '/'
    par, score_vs_par = read_validation_data(directory)
    # Calculate the best score
    best_score = calculate_best_score(par, score_vs_par)
    # Plot the validation data
    filename = directory + 'parameter_errors_formated.png'
    plot_score_vs_par(score_vs_par, par, True, filename)	
    keep_figures_live()