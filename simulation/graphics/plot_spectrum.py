'''
Plots the simulated and experimental spectra
'''

import numpy as np
import simulation.graphics.set_backend
import matplotlib.pyplot as plt
import simulation.graphics.set_style

def plot_spectrum(xs, ys, xe, ye, ranges=[], xnormalized=False, xn=[], save_figure=False, filename=''):
    fig = plt.figure(facecolor='w', edgecolor='w')
    axes = fig.gca()
    if not (xe == []):
        axes.plot(xe, ye, 'k-')
        axes.plot(xs, ys, 'r--')
        axes.legend(('exp', 'fit'), loc='upper right', frameon=False)
        if not (ranges == []):
            axes.set_xlim(-ranges['f_max'], ranges['f_max'])
            axes.set_ylim(0.0, ranges['spc_max']+0.1)
        else:
            axes.set_xlim(np.amin(xe), np.amax(xe))
            axes.set_ylim(0.0, np.amax(ye)+0.1)
        axes.set_xlabel(r'Frequency (MHz)')
        axes.set_ylabel('Amplitude')
    else:
        if xnormalized:
            axes.plot(xn, ys, 'k-')
            axes.set_xlim(np.amin(xn), np.amax(xn))
            axes.set_ylim(0, np.amax(ys)+0.1)
            axes.set_xlabel(r'$\nu_{dd}$ ($\nu_{0}$)')
            axes.set_ylabel('Amplitude')
        else:
            axes.plot(xs, ys, 'k-')
            axes.set_xlim(np.amin(xs), np.amax(xs))
            axes.set_ylim(0.0, np.amax(ys)+0.1)
            axes.set_xlabel(r'Frequency (MHz)')
            axes.set_ylabel('Amplitude')
    plt.tight_layout()
    plt.draw()
    plt.show(block=False)
    if save_figure:
        plt.savefig(filename, format='png', dpi=600)