'''
Genetic Algorithm: Plot the score in dependence of individual parameters
'''

import sys
import numpy as np
import scipy
import fitting.graphics.set_backend
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import fitting.graphics.set_style
from supplement.constants import const
import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")

# No. of plots: 1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16
alignement = [[1,1], [1,2], [1,3], [2,2], [2,3], [2,3], [2,4], [2,4], [3,3], [3,4], [3,4], [3,4], [4,4], [4,4], [4,4], [4,4]]


def plot_1d(fig, score_vs_par, var):
    # Read out the values of a gene and the corresponding score values
    name1 = var[0]
    x1 = [(x / const['variableScales'][name1]) for x in score_vs_par[name1]]
    y = score_vs_par['score']
    # Set the maximum and the minimum of y
    ymin = np.nanmin(y)
    ymax = 2 * ymin
    # Plot the figure
    axes = fig.gca()
    im = axes.scatter(x1, y, c=y, cmap='jet_r', vmin=ymin, vmax=ymax)
    axes.set_xlim(round(np.amin(x1),1), round(np.amax(x1),1))
    axes.set_xlabel(const['variableLabels'][name1])
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
    x1 = [(x / const['variableScales'][name1]) for x in score_vs_par[name1]]	
    x2 = [(x / const['variableScales'][name2]) for x in score_vs_par[name2]]	
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
    axes.set_xlabel(const['variableLabels'][name1])
    axes.set_ylabel(const['variableLabels'][name2])
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