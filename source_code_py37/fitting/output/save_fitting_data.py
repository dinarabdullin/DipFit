'''
Genetic Algorithm: Saves the results of the fitting
'''

import sys
from simulation.output.save_spectrum import save_spectrum
from simulation.output.save_timetrace import save_timetrace
from fitting.output.save_score import save_score
from fitting.output.save_fitting_parameters import save_fitting_parameters
from fitting.output.save_score_vs_par import save_score_vs_par


def save_fitting_data(fit, exp, fitPar, calc, output):
    if output['save_data']:
        # Output status
        sys.stdout.write('Saving the fitting results into the directory:\n')
        sys.stdout.write(output['directory'])
        # Save the fit to the experimental spectrum
        filename = output['directory'] + 'fit.dat'
        if calc['fitted_data'] == 'spectrum':
            save_spectrum(exp['f'], fit.best_fit, exp['f'], exp['spc'], filename)
        elif calc['fitted_data'] == 'timetrace':
            save_timetrace(exp['t'], fit.best_fit, exp['t'], exp['sig'], filename)
        # Saves the score vs optimization step
        filename = output['directory'] + 'score.dat'
        save_score(fit.best_score, filename)
        # Save the optimized fitting parameters
        filename = output['directory'] + 'parameters.dat'
        save_fitting_parameters(fit.best_parameters, filename)
        # Save the score vs individual fitting parameters
        save_score_vs_par(fit.score_vs_parameters, fitPar['errors']['variables'], output['directory'])
        # Output status
        sys.stdout.write(' [DONE]\n\n')  