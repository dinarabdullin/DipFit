'''
Genetic Algorithm: Plot the results of the fitting
'''

import sys
import numpy as np
from fitting.graphics.plot_score_vs_par import plot_score_vs_par

def plot_validation_data(fit, exp, fitPar, calc, output):
    sys.stdout.write('Plotting the validation results... ')
    # Plot the score vs genes
    filename = output['directory'] + 'parameter_errors.png'
    plot_score_vs_par(fit.score_vs_parameters, fitPar['errors']['variables'], fitPar['variables']['indices'], fitPar['variables']['bounds'], output['save_figures'], filename)	
    sys.stdout.write('[DONE]\n\n')