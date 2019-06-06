'''
Main file of the program DipFit
'''

import argparse
from input.read_config import read_config
from simulation.simulation import Simulator
from simulation.graphics.plot_simulation_data import plot_simulation_data
from simulation.output.save_simulation_data import save_simulation_data
from fitting.genetic_algorithm import GeneticAlgorithm
from fitting.graphics.plot_fitting_data import plot_fitting_data
from fitting.graphics.plot_validation_data import plot_validation_data
from fitting.output.save_fitting_data import save_fitting_data
from fitting.output.save_validation_data import save_validation_data
from supplement.make_output_directory import make_output_directory 
from supplement.keep_figures_live import keep_figures_live 


if __name__ == '__main__':
    # Read out the path to the config file
    parser = argparse.ArgumentParser()
    parser.add_argument('filepath', help="The path to the configuration file")
    args = parser.parse_args()
    configPath = args.filepath

    # Read out the config file 
    exp, spinA, spinB, mode, simPar, fitPar, calc, solution, output = read_config(configPath)

    # Make an output directory
    make_output_directory(output, configPath)

    # Spectral simulator
    sim = Simulator(calc)

    # Simulation
    if mode['simulation']:
        # Run the simulation
        sim.run_simulation(simPar, exp, spinA, spinB, calc)	
        # Save simulation results
        save_simulation_data(sim, exp, simPar['settings'], output)	
        # Plot simulation results 
        plot_simulation_data(sim, exp, simPar['settings'], calc, output)

    # Fitting
    if mode['fitting']:
        # Init the fitting mode of the simulator
        sim.init_fitting(fitPar, exp, spinA, spinB, calc)
        # Optimizer
        fit = GeneticAlgorithm(fitPar['settings'], exp, calc)
        # Run the fitting
        fit.run_optimization(fitPar, sim, exp, spinA, spinB, calc)
        # Plot fitting results
        plot_fitting_data(fit, exp, fitPar, calc, output)
        # Save fitting results
        save_fitting_data(fit, exp, fitPar, calc, output)
    
    # Validation
    if mode['validation']:
        # Init the fitting mode of the simulator
        sim.init_fitting(fitPar, exp, spinA, spinB, calc)
        # Optimizer
        fit = GeneticAlgorithm(fitPar['settings'], exp, calc)
        # Run the validation
        fit.validation(solution, fitPar, sim, exp, spinA, spinB, calc)
        # Plot validation results
        plot_validation_data(fit, exp, fitPar, calc, output)
        # Save validation results
        save_validation_data(fit, exp, fitPar, calc, output)
        
    # Keep all figures live
    keep_figures_live()