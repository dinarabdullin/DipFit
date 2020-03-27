'''
Save fitting results
'''

import sys
from fitting.output.save_fit import save_fit
from fitting.output.save_score import save_score
from fitting.output.save_fitting_parameters import save_fitting_parameters


def save_fitting_data(optimizer, expData, fitSettings, outputSettings):
    if outputSettings['save_data']:
        sys.stdout.write('Saving fitting results into the directory:\n')
        sys.stdout.write(outputSettings['directory'])
        # Save the fit
        filename = outputSettings['directory'] + 'fit.dat'
        save_fit(optimizer.best_fit, expData, fitSettings['settings']['fitted_data'], filename)
        # Saves the score vs optimization step
        filename = outputSettings['directory'] + 'score.dat'
        save_score(optimizer.best_score, filename)
        # Save the optimized fitting parameters
        filename = outputSettings['directory'] + 'parameters.dat'
        save_fitting_parameters(optimizer.best_parameters, filename)
        # Output status
        sys.stdout.write(' [DONE]\n\n')  