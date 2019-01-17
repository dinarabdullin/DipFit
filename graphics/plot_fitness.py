'''
Genetic Algorithm: Plots the fitness as a function of the number of generations
'''

import numpy as np
import set_style
import matplotlib.pyplot as plt

def plot_fitness(fitness, save_figure=False, filename=''):
	y = [v for v in fitness if not v==0]
	x = np.linspace(1,len(y),len(y))
	fig = plt.figure(facecolor='w', edgecolor='w')
	axes = fig.gca()
	axes.semilogy(x, y, linestyle='-', marker='o', color='k')
	axes.set_xlim(0, x[-1] + 1)
	plt.xlabel('The number of optimization steps')
	plt.ylabel('RMSD')	
	plt.grid(True)
	plt.tight_layout()
	plt.draw()
	plt.show(block=False)
	if save_figure:
		plt.savefig(filename, format='png', dpi=600)
	return [fig, axes]
	
def update_fitness_plot(axes, fitness):
	y = [v for v in fitness if not v==0]
	x = np.linspace(1,len(y),len(y))
	axes.clear()
	axes.semilogy(x, y, linestyle='-', marker='o', color='k')
	axes.set_xlim(0, x[-1] + 1)
	plt.xlabel('The number of optimization steps')
	plt.ylabel('RMSD')	
	plt.grid(True)
	plt.tight_layout()
	plt.draw()
	plt.show(block=False)

def close_fitness_plot(fig):
	plt.close(fig)