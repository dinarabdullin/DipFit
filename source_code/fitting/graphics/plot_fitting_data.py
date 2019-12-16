'''
Genetic Algorithm: Plot the results of the fitting
'''

import sys
import numpy as np
from fitting.graphics.plot_fit import plot_fit
from fitting.graphics.plot_score import plot_score
from fitting.graphics.plot_score_vs_par import plot_score_vs_par


def plot_fitting_data(optimizer, expData, fitSettings, valSettings, calcSettings, outputSettings):
    sys.stdout.write('Plotting the fitting results... ')
    filename = outputSettings['directory'] + 'fit.png'
    fig_fit, graph_fit = plot_fit(optimizer.best_fit, expData, fitSettings['settings']['fitted_data'], calcSettings, outputSettings['save_figures'], filename)
    
    filename = outputSettings['directory'] + 'score.png'
    fig_score, axes_score = plot_score(optimizer.best_score, outputSettings['save_figures'], filename)
    
    filename = outputSettings['directory'] + 'parameter_errors.png'
    plot_score_vs_par(optimizer.score_vs_parameters, valSettings['variables'], fitSettings, valSettings['display_threshold'], optimizer.score_threshold, outputSettings['save_figures'], filename)	
    sys.stdout.write('[DONE]\n\n')