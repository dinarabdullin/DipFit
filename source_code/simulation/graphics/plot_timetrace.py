'''
Plot the simulated and experimental time traces
'''

import numpy as np
import simulation.graphics.set_backend
import matplotlib.pyplot as plt
import simulation.graphics.set_style


def plot_timetrace(xs, ys, xe, ye, save_figure=False, filename=''):
    fig = plt.figure(facecolor='w', edgecolor='w')
    axes = fig.gca()
    if not (xe == []):    
        axes.plot(xe, ye, 'k-')
        axes.plot(xs, ys, 'r--')	
        axes.legend(('exp', 'fit'), loc='upper right', frameon=False)	
    else:
        axes.plot(xs, ys, 'r-')	
    plt.xlim([min(xs), max(xs)])
    plt.ylim([np.amin(ys)-0.1, 1.1])
    plt.xlabel(r'$\mathit{t}$ ($\mathit{\mu s}$)')
    plt.ylabel('Echo intensity (a.u.)')
    plt.tight_layout()
    plt.draw()
    plt.show(block=False)
    if save_figure:
        plt.savefig(filename, format='png', dpi=600)