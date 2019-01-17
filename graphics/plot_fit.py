'''
Genetic Algorithm: Plots the fit to the experimental spectrum
'''

import set_style
import numpy as np
import matplotlib.pyplot as plt

def plot_fit(xs, ys, xe, ye, ranges=[], save_figure=False, filename=''):      
	fig = plt.figure(facecolor='w', edgecolor='w')
	axes = fig.gca()
	axes.plot(xe, ye, 'k-')
	graph = axes.plot(xs, ys, 'r--')
	axes.legend(('exp', 'fit'), loc='upper right', frameon=False)
	if not (ranges == []):
		axes.set_xlim(-ranges['f_max'], ranges['f_max'])
		axes.set_ylim(0.0, ranges['spc_max']+0.1)
	else:
		axes.set_xlim(np.amin(xe), np.amax(xe))
		axes.set_ylim(0.0, np.amax(ye)+0.1)
	axes.set_xlabel(r'Frequency (MHz)')
	axes.set_ylabel('Amplitude')
	plt.tight_layout()
	plt.draw()
	plt.show(block=False)
	if save_figure:
		plt.savefig(filename, format='png', dpi=600)
	return [fig, graph]
		
def update_fit_plot(fig, graph, ys):
    graph[0].set_ydata(ys)
    fig.canvas.draw()
    fig.canvas.flush_events()

def close_fit_plot(fig):
	plt.close(fig)