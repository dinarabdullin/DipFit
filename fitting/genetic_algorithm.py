'''
Genetic algorithm
'''

import sys
import time
import datetime
import numpy as np
from chromosome import Chromosome
from generation import Generation
from score_function import get_fit
from get_parameters import get_parameters
from error_estimation import calculate_score_vs_parameters
from fitting.graphics.plot_fit import plot_fit, update_fit_plot, close_fit_plot
from fitting.graphics.plot_score import plot_score, update_score_plot, close_score_plot
from parameters2genes import parameters2genes
from supplement.constants import const	

class GeneticAlgorithm:

    def __init__(self, settings, exp, calc):
        self.num_generations = settings['num_generations']
        self.size_generation = settings['size_generation']
        self.prob_crossover = settings['prob_crossover']
        self.prob_mutation = settings['prob_mutation']
        self.display_graphics = settings['display_graphics']
        if calc['fitted_data'] == 'spectrum':
            self.best_fit = np.zeros(exp['spc'].size)  
        elif calc['fitted_data'] == 'timetrace':
            self.best_fit = np.zeros(exp['sig'].size)  
        self.best_score = np.zeros(self.num_generations)
        self.best_parameters = {}
        self.score_vs_parameters = []

    def run_optimization(self, fitPar, sim, exp, spinA, spinB, calc):
        # Display status
        sys.stdout.write('Starting fitting...\n')
        time_start = time.time()  
        # Increment over generations
        for i in range(self.num_generations):     
            if (i == 0):
                # Create the first generation
                generation = Generation(self.size_generation)
                generation.first_generation(fitPar['variables']['bounds'])
            else:
                # Create the next generation
                generation.produce_offspring(fitPar['variables']['bounds'], self.prob_crossover, self.prob_mutation)
            # Score the generation
            generation.score_chromosomes(fitPar['variables']['indices'], fitPar['variables']['fixed'], sim, exp, spinA, spinB, calc)
            # Sort chromosomes according to their score
            generation.sort_chromosomes() 
            # Save the best score in each optimization step
            self.best_score[i] = generation.chromosomes[0].score
            # Display some graphics
            if self.display_graphics:
                # Calculate the best fit so far
                best_fit = get_fit(generation.chromosomes[0].genes, fitPar['variables']['indices'], fitPar['variables']['fixed'], sim, exp, spinA, spinB, calc)
                # Plot the experimental and simulated spectra
                if (i == 0):
                    if calc['fitted_data'] == 'spectrum':
                        # Plot the dipolar spectrum
                        fig_fit, graph_fit = plot_fit(exp['f'], best_fit, exp['f'], exp['spc'], calc)    
                    elif calc['fitted_data'] == 'timetrace':
                        # Plot the dipolar timetrace
                        fig_fit, graph_fit = plot_fit(exp['t'], best_fit, exp['t'], exp['sig'], calc)    
                    fig_score, axes_score = plot_score(self.best_score)
                elif ((i > 0) and (i < (self.num_generations - 1))):
                    update_fit_plot(fig_fit, graph_fit, best_fit)
                    update_score_plot(axes_score, self.best_score)
                elif (i == (self.num_generations-1)):
                    close_fit_plot(fig_fit)
                    close_score_plot(fig_score)
            # Output the current status
            sys.stdout.write('\r')
            sys.stdout.write("Optimization step %d / %d: RMSD = %f" % (i+1, self.num_generations, self.best_score[i]))
            sys.stdout.flush()
        # Calculate the best fit
        if self.display_graphics:
            self.best_fit = best_fit
        else:	
            self.best_fit = get_fit(generation.chromosomes[0].genes, fitPar['variables']['indices'], fitPar['variables']['fixed'], sim, exp, spinA, spinB, calc)
        # Store the best genes
        self.best_parameters = get_parameters(generation.chromosomes[0].genes, fitPar['variables']['indices'], fitPar['variables']['fixed'])
        # Display status
        time_finish = time.time()
        time_elapsed = str(datetime.timedelta(seconds = time_finish - time_start))
        sys.stdout.write('\n')
        sys.stdout.write('Optimization is finished. Total duration: %s' % (time_elapsed))
        sys.stdout.write('\n\n')	
        # Estimate the precision of the genes obtained
        if not (fitPar['errors']['variables'] == []):	
            sys.stdout.write('Estimating the precision of fitting parameters... \n')
            time_start = time.time()
            self.score_vs_parameters, parameter_errors = calculate_score_vs_parameters(self.best_parameters, fitPar, sim, exp, spinA, spinB, calc)
            self.best_parameters = get_parameters(generation.chromosomes[0].genes, fitPar['variables']['indices'], fitPar['variables']['fixed'], parameter_errors)
            time_finish = time.time()
            time_elapsed = str(datetime.timedelta(seconds = time_finish - time_start))
            sys.stdout.write('Error estimation is finished. Total duration: %s' % (time_elapsed))
            sys.stdout.write('\n\n')
        # Make a summary of the fitting parameters
        sys.stdout.write('Optimized parameters:\n')
        sys.stdout.write("{0:<16s} {1:<16s} {2:<16s} {3:<16s}\n".format('Parameter', 'Value', 'Optimized', 'Precision (+/-)'))
        for name in const['variableNames']:
            parameter = self.best_parameters[name]
            sys.stdout.write("{0:<16s} {1:<16.3f} {2:<16s} {3:<16.3f}\n".format(parameter['longname'], parameter['value']/const['variableScales'][name], parameter['optimized'], parameter['precision']/const['variableScales'][name]))
        sys.stdout.write('\n')

    def validation(self, solution, fitPar, sim, exp, spinA, spinB, calc):
        # Display status
        sys.stdout.write('Starting validation...\n')
        time_start = time.time()  
        # Store the best genes
        self.best_parameters = solution
        # Estimate the precision of the genes obtained
        if not (fitPar['errors']['variables'] == []):	
            self.score_vs_parameters, parameter_errors = calculate_score_vs_parameters(self.best_parameters, fitPar, sim, exp, spinA, spinB, calc)
            best_genes = parameters2genes(self.best_parameters)
            self.best_parameters = get_parameters(best_genes, fitPar['variables']['indices'], fitPar['variables']['fixed'], parameter_errors)
        # Display status
        time_finish = time.time()
        time_elapsed = str(datetime.timedelta(seconds = time_finish - time_start))
        sys.stdout.write('\n')
        sys.stdout.write('Validation is finished. Total duration: %s' % (time_elapsed))
        sys.stdout.write('\n\n')
        # Make a summary of the fitting parameters
        sys.stdout.write('Optimized parameters:\n')
        sys.stdout.write("{0:<16s} {1:<16s} {2:<16s} {3:<16s}\n".format('Parameter', 'Value', 'Optimized', 'Precision (+/-)'))
        for name in const['variableNames']:
            parameter = self.best_parameters[name]
            sys.stdout.write("{0:<16s} {1:<16.3f} {2:<16s} {3:<16.3f}\n".format(parameter['longname'], parameter['value']/const['variableScales'][name], parameter['optimized'], parameter['precision']/const['variableScales'][name]))
        sys.stdout.write('\n')
        