'''
Genetic Algorithm: Plots the score as a function of individual genes
'''

import fitting.graphics.set_style

import sys
import numpy as np
import scipy
import fitting.graphics.set_backend
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
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

def plot_1d(fig, data, var, indices=[], bounds=[]):
    # Read out the values of a gene and the corresponding score values
    name1 = var[0]
    x1 = [(x / const['variableScales'][name1]) for x in data[name1]]
    y = data['score']
    # Set the limits for y
    ymin = np.min(y)
    ymax = 2.0 * ymin
    # Plot the graph
    axes = fig.gca()
    im = axes.scatter(x1, y, c=y, cmap='jet_r', vmin=ymin, vmax=ymax, linewidth='0')
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

def plot_2d(fig, data, var, indices=[], bounds=[]):
    # Read out the values of a pair of genes and the corresponding score values
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
    X, Y = np.meshgrid(x1r, x2r)
    Z = griddata(x1, x2, y, x1r, x2r, interp='linear')
    # Set the limits for y
    ymin = np.min(Z)
    ymax = 2.0 * ymin
    # Plot the graph
    axes = fig.gca()
    im = axes.pcolor(X, Y, Z, cmap='jet_r', vmin=ymin, vmax=ymax)
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
    return im

def plot_score_vs_par(score_vs_par, par, indices=[], bounds=[], save_figure=False, filename=''): 
	M = len(par)
	Mp = count_plots(par)
	Np = 1
	fig = plt.figure(figsize=(18,9), facecolor='w', edgecolor='w')
	for i in range(M):
		dim = len(par[i])
		if (dim == 1):
			plt.subplot(alignement[Mp-1][0], alignement[Mp-1][1], Np)
			Np = Np + 1
			im = plot_1d(fig, score_vs_par[i], par[i], indices, bounds)
		elif (dim == 2):
			plt.subplot(alignement[Mp-1][0], alignement[Mp-1][1], Np)
			Np = Np + 1
			im = plot_2d(fig, score_vs_par[i], par[i], indices, bounds)
	plt.tight_layout()
	plt.subplots_adjust(bottom=0.1, right=0.85, top=0.9)
	cax = plt.axes([0.9, 0.3, 0.02, 0.4]) # left, bottom, width, height
	plt.colorbar(im, cax=cax, orientation='vertical')
	plt.text(0.8, 1.1, 'RMSD')
	plt.draw()
	plt.show(block=False)
	if save_figure:
		plt.savefig(filename, format='png', dpi=600)