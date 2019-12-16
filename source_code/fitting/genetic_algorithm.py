'''
Genetic algorithm
'''

import sys
import time
import datetime
import numpy as np
from fitting.chromosome import Chromosome
from fitting.generation import Generation
from fitting.score_function import get_fit
from fitting.get_parameters import get_parameters
from fitting.error_estimation import calculate_score_vs_parameters, calculation_error
from fitting.graphics.plot_fit import plot_fit, update_fit_plot, close_fit_plot
from fitting.graphics.plot_score import plot_score, update_score_plot, close_score_plot
from fitting.parameters2genes import parameters2genes
from supplement.constants import const	


class GeneticAlgorithm:
    def __init__(self, settings, expData):
        self.num_generations = settings['num_generations']
        self.size_generation = settings['size_generation']
        self.prob_crossover = settings['prob_crossover']
        self.prob_mutation = settings['prob_mutation']
        self.fitted_data = settings['fitted_data']
        self.display_graphics = settings['display_graphics']
        if self.fitted_data == 'spectrum':
            self.best_fit = np.zeros(expData['spc'].size)  
        elif self.fitted_data == 'timetrace':
            self.best_fit = np.zeros(expData['sig'].size)  
        self.best_score = np.zeros(self.num_generations)
        self.best_parameters = {}
        self.score_vs_parameters = []
        self.calc_error = 0
        self.score_threshold = 0

    def run_optimization(self, fitSettings, simulator, expData, spinA, spinB, calcSettings):
        sys.stdout.write('Starting fitting...\n')
        time_start = time.time()  
        for i in range(self.num_generations):     
            if (i == 0):
                # Create the first generation
                generation = Generation(self.size_generation)
                generation.first_generation(fitSettings['variables']['bounds'])
            else:
                # Create the next generation
                generation.produce_offspring(fitSettings['variables']['bounds'], self.prob_crossover, self.prob_mutation)
            # Score the generation
            generation.score_chromosomes(fitSettings, simulator, expData, spinA, spinB, calcSettings)
            # Sort chromosomes according to their score
            generation.sort_chromosomes() 
            # Save the best score in each optimization step
            self.best_score[i] = generation.chromosomes[0].score
            # Display some graphics
            if self.display_graphics:
                # Calculate the best fit so far
                best_fit = get_fit(generation.chromosomes[0].genes, fitSettings, simulator, expData, spinA, spinB, calcSettings)
                if (i == 0):
                    fig_fit, graph_fit = plot_fit(best_fit, expData, self.fitted_data, calcSettings)   
                    fig_score, axes_score = plot_score(self.best_score)
                elif ((i > 0) and (i < self.num_generations-1)):
                    update_fit_plot(fig_fit, graph_fit, best_fit)
                    update_score_plot(axes_score, self.best_score)
                elif (i == self.num_generations-1):
                    close_fit_plot(fig_fit)
                    close_score_plot(fig_score)                
            sys.stdout.write('\r')
            sys.stdout.write("Optimization step %d / %d: RMSD = %f" % (i+1, self.num_generations, self.best_score[i]))
            sys.stdout.flush()
        # Calculate the best fit
        if self.display_graphics:
            self.best_fit = best_fit
        else:	
            self.best_fit = get_fit(generation.chromosomes[0].genes, fitSettings, simulator, expData, spinA, spinB, calcSettings)
        # Store the best genes
        self.best_parameters = get_parameters(generation.chromosomes[0].genes, fitSettings['variables']['indices'], fitSettings['variables']['fixed'])
        time_finish = time.time()
        time_elapsed = str(datetime.timedelta(seconds = time_finish - time_start))
        sys.stdout.write('\n')
        sys.stdout.write('Fitting is finished. Total duration: %s' % (time_elapsed))
        sys.stdout.write('\n\n')	
         
    def validation(self, valSettings, fitSettings, simulator, expData, spinA, spinB, calcSettings):   
        if not (valSettings['variables'] == []):
            sys.stdout.write('Validating the optimized fitting parameters... \n')
            time_start = time.time()
            # Determine the calculation error
            self.calc_error, score_min = calculation_error(self.best_parameters, fitSettings, simulator, expData, spinA, spinB, calcSettings)
            sys.stdout.write('RMSD calculation error = %f\n' % (self.calc_error))
            # Set the score threshold
            if valSettings['threshold']:
                self.score_threshold = valSettings['threshold'] * score_min
                score_error = (valSettings['threshold'] - 1.0) * score_min
                if (score_error < self.calc_error):
                    sys.stdout.write('Warning: The error threshold is below the calculation error!\n')
                    self.score_threshold = score_min + self.calc_error
            else:
                self.score_threshold = score_min + self.calc_error
            self.score_vs_parameters, parameter_errors = calculate_score_vs_parameters(self.best_parameters, valSettings, fitSettings, simulator, expData, spinA, spinB, calcSettings, self.score_threshold)
            best_genes = parameters2genes(self.best_parameters)
            self.best_parameters = get_parameters(best_genes, fitSettings['variables']['indices'], fitSettings['variables']['fixed'], parameter_errors)
            time_finish = time.time()
            time_elapsed = str(datetime.timedelta(seconds = time_finish - time_start))
            sys.stdout.write('Validation is finished. Total duration: %s' % (time_elapsed))
            sys.stdout.write('\n\n')
            
    def print_optimized_parameters(self):
        sys.stdout.write('Optimized fitting parameters:\n')
        sys.stdout.write("{0:<16s} {1:<16s} {2:<16s} {3:<16s}\n".format('Parameter', 'Value', 'Optimized', 'Precision (+/-)'))
        for name in const['variableNames']:
            parameter = self.best_parameters[name]
            sys.stdout.write("{0:<16s} {1:<16.3f} {2:<16s} {3:<16.3f}\n".format(parameter['longname'], parameter['value']/const['variableScales'][name], parameter['optimized'], parameter['precision']/const['variableScales'][name]))
        sys.stdout.write('\n')
        