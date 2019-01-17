'''
Genetic Algorithm: Plots the fitness as a function of individual genes
'''

import numpy as np
import set_style
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import sys
import warnings
if not sys.warnoptions:
    warnings.simplefilter("ignore")
from constants import const

# No. of plots: 1      2      3      4      5      6      7      8      9     10     11     12     13     14     15     16
alignement = [[1,1], [1,2], [1,3], [2,2], [2,3], [2,3], [2,4], [2,4], [3,3], [3,4], [3,4], [3,4], [4,4], [4,4], [4,4], [4,4]]

def plot_fitness_vs_genes(fitness_vs_genes, gene_sets, indices=[], bounds=[], save_figure=False, filename=''):
	M = len(gene_sets)
	fig = plt.figure(figsize=(18,9), facecolor='w', edgecolor='w')
	for i in range(M):
		plt.subplot(alignement[M-1][0], alignement[M-1][1], i+1)
		dim = len(gene_sets[i])
		if (dim == 1):
			# Read out the values of a gene and the corresponding fitness values
			name1 = gene_sets[i][0]
			x1 = [(x / const['variableScales'][name1]) for x in fitness_vs_genes[i][name1]]			
			y = fitness_vs_genes[i]['fitness']
			# Set the limits for y
			ymin = np.min(y)
			ymax = 2.0 * ymin
			# Plot the graph
			axes = fig.gca()
			axes.scatter(x1, y, c=y, cmap='jet_r', vmin=ymin, vmax=ymax, linewidth='0')
			if not (indices == []) and (bounds == []):
				axes.set_xlim([bounds[indices[name1]][0]/const['variableScales'][name1], bounds[indices[name1]][1]/const['variableScales'][name1]])
			else:
				axes.set_xlim(np.amin(x1), np.amax(x1))
			axes.set_xlabel(const['variableLabels'][name1])
			axes.set_ylabel('RMSD')
		elif (dim == 2):
			# Read out the values of a pair of genes and the corresponding fitness values
			name1 = gene_sets[i][0]
			name2 = gene_sets[i][1]
			x1 = [(x / const['variableScales'][name1]) for x in fitness_vs_genes[i][name1]]	
			x2 = [(x / const['variableScales'][name2]) for x in fitness_vs_genes[i][name2]]	
			y = fitness_vs_genes[i]['fitness']
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
		x0,x1 = axes.get_xlim()
		y0,y1 = axes.get_ylim()
		axes.set_aspect((x1-x0)/(y1-y0))
	plt.tight_layout()
	plt.subplots_adjust(bottom=0.1, right=0.85, top=0.9)
	cax = plt.axes([0.9, 0.3, 0.02, 0.4]) # left, bottom, width, height
	plt.colorbar(im, cax=cax, orientation='vertical')
	plt.text(0.8, 1.1, 'RMSD')
	plt.draw()
	plt.show(block=False)
	if save_figure:
		plt.savefig(filename, format='png', dpi=600)