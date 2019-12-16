'''
Genetic Algorithm: Plot the results of the fitting
'''

import sys
import numpy as np
from fitting.graphics.plot_score_vs_par import plot_score_vs_par


def plot_validation_data(optimizer, valSettings, fitSettings, outputSettings):
    sys.stdout.write('Plotting the validation results... ')
    filename = outputSettings['directory'] + 'parameter_errors.png'
    plot_score_vs_par(optimizer.score_vs_parameters, valSettings['variables'], fitSettings, valSettings['display_threshold'], optimizer.score_threshold, outputSettings['save_figures'], filename)	
    sys.stdout.write('[DONE]\n\n')