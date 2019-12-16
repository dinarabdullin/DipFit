'''
Saves the validation results
'''

import sys
from fitting.output.save_fitting_parameters import save_fitting_parameters
from fitting.output.save_score_vs_par import save_score_vs_par


def save_validation_data(optimizer, valSettings, outputSettings):
    if outputSettings['save_data']:
        sys.stdout.write('Saving the fitting results into the directory:\n')
        sys.stdout.write(outputSettings['directory'])
        # Save the optimized fitting parameters
        filename = outputSettings['directory'] + 'parameters.dat'
        save_fitting_parameters(optimizer.best_parameters, filename)
        # Save the score vs individual fitting parameters
        save_score_vs_par(optimizer.score_vs_parameters, valSettings['variables'], outputSettings['directory'])
        # Output status
        sys.stdout.write(' [DONE]\n\n')  