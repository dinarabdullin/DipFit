'''
Genetic Algorithm: Plot the score in dependence of optimization step
'''

import numpy as np
import fitting.graphics.set_backend
import matplotlib.pyplot as plt
import fitting.graphics.set_style


def plot_score(score, save_figure=False, filename=''): 
	y = [v for v in score if not v==0]
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


def update_score_plot(axes, score):
	y = [v for v in score if not v==0]
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


def close_score_plot(fig):
	plt.close(fig)