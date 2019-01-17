'''
Plots the flip probability vs g-factor and temperature
'''

import numpy as np
import set_style
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def plot_flip_probabilities(g, l, T, save_figure=False, filename=''): 
	fig = plt.figure(facecolor='w', edgecolor='w')
	axes = fig.gca()
	Nt = len(T)
	colors = cm.rainbow(np.linspace(0, 1, Nt))
	for i in range(Nt):
		axes.plot(g, l[i], color=colors[i])
	axes.set_xlim(np.amin(g), np.amax(g))
	axes.set_ylim(np.amin(l)-0.1, np.amax(l)+0.1)
	axes.set_xlabel(r'$\mathit{g^\prime_{1eff}}$')
	axes.set_ylabel(r'$p$')
	labels = [(str(t) + ' K') for t in T]
	axes.legend(labels=labels, loc='lower left')
	plt.tight_layout()
	plt.draw()
	plt.show(block=False)
	if save_figure:
		plt.savefig(filename, format='png', dpi=600)