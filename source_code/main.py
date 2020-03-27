'''
Main file of DipFit
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
    # Read out the config file 
    parser = argparse.ArgumentParser()
    parser.add_argument('filepath', help="A path to a configuration file")
    args = parser.parse_args()
    configPath = args.filepath
    mode, expData, spinA, spinB, simSettings, fitSettings, valSettings, calcSettings, outputSettings = read_config(configPath)

    # Make an output directory
    make_output_directory(outputSettings, configPath)

    # Spectral simulator
    simulator = Simulator(calcSettings)

    # Simulation
    if mode['simulation']:
        # Run the simulation
        simulator.run_simulation(simSettings, expData, spinA, spinB, calcSettings)	
        # Save simulation results
        save_simulation_data(simulator, simSettings, expData, outputSettings)	
        # Plot simulation results 
        plot_simulation_data(simulator, simSettings, expData, calcSettings, outputSettings)

    # Fitting
    if mode['fitting']:
        # Init the fitting mode of the simulator
        simulator.init_fitting(fitSettings, expData, spinA, spinB, calcSettings)
        if fitSettings['settings']['method'] == "genetic":
            # Optimizer
            optimizer = GeneticAlgorithm(fitSettings['settings'], expData)
            # Run the fitting
            optimizer.run_optimization(fitSettings, simulator, expData, spinA, spinB, calcSettings)
        # Save fitting results
        save_fitting_data(optimizer, expData, fitSettings, outputSettings)
        # Plot fitting results
        plot_fitting_data(optimizer, expData, fitSettings, calcSettings, outputSettings)
        # Validate the fitting parameters
        optimizer.validation(valSettings, fitSettings, simulator, expData, spinA, spinB, calcSettings)
        # Save validation results
        save_validation_data(optimizer, valSettings, outputSettings)
        # Plot validation results
        plot_validation_data(optimizer, valSettings, fitSettings, outputSettings)
        # Display the optmized fitting parameters
        optimizer.print_optimized_parameters()    

    # Validation
    if mode['validation']:
        # Init the fitting mode of the simulator
        simulator.init_fitting(fitSettings, expData, spinA, spinB, calcSettings)
        # Optimizer
        if fitSettings['settings']['method'] == "genetic":
            optimizer = GeneticAlgorithm(fitSettings['settings'], expData)
        # Run the validation
        optimizer.best_parameters = valSettings['optimized_parameters']
        optimizer.validation(valSettings, fitSettings, simulator, expData, spinA, spinB, calcSettings)
        # Save validation results
        save_validation_data(optimizer, valSettings, outputSettings)
        # Plot validation results
        plot_validation_data(optimizer, valSettings, fitSettings, outputSettings)
        # Display the optmized fitting parameters
        optimizer.print_optimized_parameters()
           
    # Keep all figures live
    keep_figures_live()