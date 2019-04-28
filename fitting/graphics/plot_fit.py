'''
Genetic Algorithm: Plots the fit to the experimental spectrum
'''

import numpy as np
import fitting.graphics.set_backend
import matplotlib.pyplot as plt
import fitting.graphics.set_style

def plot_fit(xs, ys, xe, ye, calc, save_figure=False, filename=''):      
	fig = plt.figure(facecolor='w', edgecolor='w')
	axes = fig.gca()
	axes.plot(xe, ye, 'k-')
	graph = axes.plot(xs, ys, 'r--')
	axes.legend(('exp', 'fit'), loc='upper right', frameon=False)
	if (calc['fitted_data'] == 'spectrum'):
		if not (calc['f_max'] == 0):
			axes.set_xlim(-calc['f_max'], calc['f_max'])
		else:
			axes.set_xlim(np.amin(xe), np.amax(xe))
		axes.set_ylim(0.0, calc['spc_max']+0.1)
		axes.set_xlabel(r'Frequency (MHz)')
		axes.set_ylabel('Amplitude')
	elif (calc['fitted_data'] == 'timetrace'):
		if not (calc['t_max'] == 0):
			axes.set_xlim(calc['t_min'], calc['t_max'])
		else:
			axes.set_xlim(np.amin(xe), np.amax(xe))
		axes.set_ylim(np.amin(ye)-0.2, np.amax(ye)+0.1)
		axes.set_xlabel(r'$\mathit{t}$ ($\mathit{\mu s}$)')
		axes.set_ylabel('Echo intensity (a.u.)')
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