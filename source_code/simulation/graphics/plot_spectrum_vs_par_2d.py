'''
Plot the dependence of the simulated spectra on the parameter 'par'
'''

import numpy as np
import math
import simulation.graphics.set_backend
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import simulation.graphics.set_style


def plot_spectrum_vs_par_2d(xs, par, ys, ranges=[], xnormalized=False, xn=[], save_figure=False, filename='', par_label='', invert_paxis=False):
    fig = plt.figure(figsize=(6,10), facecolor='w', edgecolor='w')
    fig.subplots_adjust(left=.12, right=.85, bottom=.10, top=.95)
    axes = fig.gca()
    Nx = xs.size
    Npar = par.size
    Nskip = Nskip = int((Npar-1) / 10) + 1
    count = 0
    # Set the value of y-shift
    if (invert_paxis):
        yshift = -0.5 * np.amax(ys)
        baseline = (-yshift * Npar) * np.ones(Nx)
    else:
        yshift = 0.7 * np.amax(ys)
        baseline = np.zeros(Nx)
    if xnormalized:
        for i in range(Npar):
            axes.plot(xn, ys[i]+baseline, 'k-')          
            if (count == 0):
                plt.text(np.amax(xn) + 0.5, baseline[0] - 0.07, str(int(par[i])))
            count = count + 1
            if (count == Nskip):
                count = 0
            baseline = baseline + yshift * np.ones(Nx)
            if (i == Npar-1):
                plt.text(np.amax(xn) + 0.5, baseline[0], par_label)
        axes.set_xlim(np.amin(xn), np.amax(xn))
        axes.set_xlabel(r'$\nu_{dd}$ ($\nu_{0}$)')
        axes.set_xticks(np.linspace(-math.floor(np.amax(xn)), math.floor(np.amax(xn)), 2 * int(math.floor(np.amax(xn))) + 1))
        axes.xaxis.grid(color='darkgray', linestyle='--', linewidth=1)
    else:
        for i in range(Npar):
            axes.plot(xs, ys[i]+baseline, 'k-')
            if (count == 0):
                plt.text(np.amax(xs) + 0.5, baseline[0] - 0.07, str(int(par[i])))
            count = count + 1
            if (count == Nskip):
                count = 0
            baseline = baseline + yshift * np.ones(Nx)
            if (i == Npar-1):   
                plt.text(np.amax(xs) + 0.5, baseline[0] + 0.4*Nskip*yshift, par_label)
        axes.set_xlim(np.amin(xs), np.amax(xs))
        axes.set_xlabel(r'Frequency (MHz)')
    #axes.set_ylabel('Amplitude')   
    axes.spines['left'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.spines['top'].set_visible(False)
    axes.tick_params(axis = 'y', which = 'both', left = False, right = False, labelleft = False, labelright = False)
    axes.tick_params(axis = 'x', which = 'both', top = False, bottom = True)
    #plt.tight_layout()
    plt.draw()
    plt.show(block=False)
    if save_figure:
        plt.savefig(filename, format='png', dpi=600)