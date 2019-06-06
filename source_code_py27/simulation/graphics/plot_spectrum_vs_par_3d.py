'''
Plot the dependence of the simulated spectra on the parameter 'par'
'''

import numpy as np
import simulation.graphics.set_backend
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import simulation.graphics.set_style


def plot_spectrum_vs_par_3d(xs, par, ys, xnormalized=False, xn=[], save_figure=False, filename='', par_label='', invert_paxis=False): 
    # Colormap
    Npar = par.size
    Nc = int(1.5 * Npar)
    cmap = plt.cm.get_cmap('gray', Nc)
    cmaplist = [cmap(i) for i in range(Nc)]
    # Plot the spectrum vs par
    fig = plt.figure(figsize=(9,9), facecolor='w', edgecolor='w')
    axes = fig.gca(projection='3d') 
    if xnormalized:
        for i in range(Npar):
            axes.plot(xn, par[i]*np.ones(xs.size), ys[i], color=cmaplist[i])
        axes.set_xlim(np.amin(xn), np.amax(xn))
        axes.set_ylim(np.amin(par), np.amax(par))
        axes.set_zlim(0.0, np.amax(ys)+0.1)
        axes.set_xlabel(r'$\nu_{dd}$ ($\nu_{0}$)')
        axes.set_ylabel(par_label)
        axes.set_zlabel('Amplitude')
    else:
        for i in range(Npar):
            axes.plot(xs, par[i]*np.ones(xs.size), ys[i], color=cmaplist[i])
        axes.set_xlim(np.amin(xs), np.amax(xs))
        axes.set_ylim(np.amin(par), np.amax(par))
        axes.set_zlim(0.0, np.amax(ys)+0.1)
        axes.set_xlabel(r'Frequency (MHz)')
        axes.set_ylabel(par_label)
        axes.set_zlabel('Amplitude')
    axes.view_init(elev=54, azim=-85)
    if (invert_paxis):
        axes.invert_yaxis()
    plt.tight_layout()
    plt.draw()
    plt.show(block=False)
    if save_figure:
        plt.savefig(filename, format='png', dpi=600)