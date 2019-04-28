'''
Genetic Algorithm: Plot the results of the fitting
'''

import sys
import numpy as np
from fitting.graphics.plot_fit import plot_fit
from fitting.graphics.plot_score import plot_score
from fitting.graphics.plot_score_vs_par import plot_score_vs_par

def plot_fitting_data(fit, exp, fitPar, calc, output):
    sys.stdout.write('Plotting the fitting results... ')
    # Plot the fit
    filename = output['directory'] + 'fit.png'
    if calc['fitted_data'] == 'spectrum':
        fig_fit, graph_fit = plot_fit(exp['f'], fit.best_fit, exp['f'], exp['spc'], calc, output['save_figures'], filename)
    elif calc['fitted_data'] == 'timetrace':
        fig_fit, graph_fit = plot_fit(exp['t'], fit.best_fit, exp['t'], exp['sig'], calc, output['save_figures'], filename)
    # Plot the score vs the number of the generation
    filename = output['directory'] + 'score.png'
    fig_score, axes_score = plot_score(fit.best_score, output['save_figures'], filename)
    # Plot the score vs genes
    filename = output['directory'] + 'parameter_errors.png'
    plot_score_vs_par(fit.score_vs_parameters, fitPar['errors']['variables'], fitPar['variables']['indices'], fitPar['variables']['bounds'], output['save_figures'], filename)	
    sys.stdout.write('[DONE]\n\n')