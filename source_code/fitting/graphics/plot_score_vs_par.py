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


def count_plots(sets):
    Mp = 0
    for subset in sets: 
        dim = len(subset)
        if (dim == 1) or (dim == 2):
            Mp = Mp + 1
    return Mp  


def set_limits(v, use_threshold, threshold):
    if use_threshold:
        if threshold:
            vmin = threshold
            vmax = 2 * threshold
        else: 
            vmin = np.nanmin(v)
            vmax = 2 * vmin
    else:
        if threshold:
            vmin = np.nanmin(v)
            vmax = 2 * threshold
        else: 
            vmin = np.nanmin(v)
            vmax = 2 * vmin 
    return [vmin, vmax]


def plot_1d(fig, data, var, display_threshold, score_threshold, indices=[], bounds=[]):
    # Read out the values of a gene and the corresponding score values
    name1 = var[0]
    x1 = [(x / const['variableScales'][name1]) for x in data[name1]]
    y = data['score']
    # Set the limits for y
    ymin, ymax = set_limits(y, display_threshold, score_threshold)
    # Plot the graph
    axes = fig.gca()
    im = axes.scatter(x1, y, c=y, cmap='jet_r', vmin=ymin, vmax=ymax)
    if not (indices == []) and (bounds == []):
        axes.set_xlim([bounds[indices[name1]][0]/const['variableScales'][name1], bounds[indices[name1]][1]/const['variableScales'][name1]])
    else:
        axes.set_xlim(np.amin(x1), np.amax(x1))
    axes.set_xlabel(const['variableLabels'][name1])
    axes.set_ylabel('RMSD')
    # Make the axis of equal scale
    xl, xh = axes.get_xlim()
    yl, yh = axes.get_ylim()
    axes.set_aspect( (xh-xl)/(yh-yl) )
    return im


def plot_2d(fig, data, var, display_threshold, score_threshold, indices=[], bounds=[]):
    # Read out the values of variables
    name1 = var[0]
    name2 = var[1]
    x1 = [(x / const['variableScales'][name1]) for x in data[name1]]	
    x2 = [(x / const['variableScales'][name2]) for x in data[name2]]	
    y = data['score']
    # Interpolate the data on a regular grid
    x1min = np.min(x1)
    x1max = np.max(x1)
    x2min = np.min(x2)
    x2max = np.max(x2)
    x1r = np.linspace(x1min, x1max, num=200) 
    x2r = np.linspace(x2min, x2max, num=200)
    X, Y = np.mgrid[x1min:x1max:200j, x2min:x2max:200j]
    Z = griddata((x1, x2), y, (X, Y), method='linear')
    # Set the limits for Z
    zmin, zmax = set_limits(Z, display_threshold, score_threshold)
    # Plot the figure
    axes = fig.gca()
    im = axes.pcolor(X, Y, Z, cmap='jet_r', vmin=zmin, vmax=zmax)
    if not (indices == []) and (bounds == []):
        axes.set_xlim([bounds[indices[name1]][0]/const['variableScales'][name1], bounds[indices[name1]][1]/const['variableScales'][name1]])
        axes.set_ylim([bounds[indices[name2]][0]/const['variableScales'][name2], bounds[indices[name2]][1]/const['variableScales'][name2]])
    else:
        axes.set_xlim(np.amin(x1), np.amax(x1))
        axes.set_ylim(np.amin(x2), np.amax(x2))
    axes.set_xlabel(const['variableLabels'][name1])
    axes.set_ylabel(const['variableLabels'][name2])
    # Make the axis of equal scale
    xl, xh = axes.get_xlim()
    yl, yh = axes.get_ylim()
    axes.set_aspect( (xh-xl)/(yh-yl) )
    # Make ticks
    plt.xticks(np.linspace(bounds[indices[name1]][0]/const['variableScales'][name1], bounds[indices[name1]][1]/const['variableScales'][name1], 3))
    plt.yticks(np.linspace(bounds[indices[name2]][0]/const['variableScales'][name2], bounds[indices[name2]][1]/const['variableScales'][name2], 3))
    axes.xaxis.labelpad = 0
    axes.yaxis.labelpad = 0
    return im


def plot_score_vs_par(score_vs_par, par, fitSettings, display_threshold, score_threshold, save_figure=False, filename=''): 
    indices = fitSettings['variables']['indices']
    bounds = fitSettings['variables']['bounds'] 
    M = len(par)
    Mp = count_plots(par)
    Np = 1
    fig = plt.figure(figsize=(18,9), facecolor='w', edgecolor='w')
    for i in range(M):
        dim = len(par[i])
        if (dim == 1):
            plt.subplot(alignement[Mp-1][0], alignement[Mp-1][1], Np)
            Np = Np + 1
            im = plot_1d(fig, score_vs_par[i], par[i], display_threshold, score_threshold, indices, bounds)
        elif (dim == 2):
            plt.subplot(alignement[Mp-1][0], alignement[Mp-1][1], Np)
            Np = Np + 1
            im = plot_2d(fig, score_vs_par[i], par[i], display_threshold, score_threshold, indices, bounds)
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15, right=0.80, top=0.9)
    cax = plt.axes([0.85, 0.3, 0.02, 0.4]) # left, bottom, width, height  
    plt.colorbar(im, cax=cax, orientation='vertical')
    plt.text(0.90, 1.05, 'RMSD', transform=cax.transAxes)
    plt.draw()
    plt.show(block=False)
    if save_figure:
        plt.savefig(filename, format='png', dpi=600)