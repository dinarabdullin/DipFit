'''
Genetic Algorithm: Saves the results of the fitting
'''

import sys
from fitting.output.save_fitting_parameters import save_fitting_parameters
from fitting.output.save_score_vs_par import save_score_vs_par


def save_validation_data(fit, exp, fitPar, calc, output):
    if output['save_data']:
        # Output status
        sys.stdout.write('Saving the fitting results into the directory:\n')
        sys.stdout.write(output['directory'])
        # Save the optimized fitting parameters
        filename = output['directory'] + 'parameters.dat'
        save_fitting_parameters(fit.best_parameters, filename)
        # Save the score vs individual fitting parameters
        save_score_vs_par(fit.score_vs_parameters, fitPar['errors']['variables'], output['directory'])
        # Output status
        sys.stdout.write(' [DONE]\n\n')  