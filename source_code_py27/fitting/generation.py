'''
Generation class
'''

import sys
from copy import deepcopy
import numpy as np
from multiprocessing import Pool
from functools import partial
from chromosome import Chromosome
from score_function import score_function


class Generation:
	def __init__(self, size_generation):
		self.size = size_generation
		self.chromosomes = []
		
	def first_generation(self, bounds):
		for i in range(self.size):
			chromosome = Chromosome(bounds)
			self.chromosomes.append(chromosome)

	def tournament_selection(self):
		idx1 = np.random.random_integers(low=0, high=self.size - 1)
		idx2 = np.random.random_integers(low=0, high=self.size - 1)
		if self.chromosomes[idx1].score < self.chromosomes[idx2].score:
			return idx1
		else:
			return idx2

	def crossover_chromosomes(self, chromosome1, chromosome2, prob_crossover):
		if (np.random.rand() <= prob_crossover):
			# Store the genes of chromosome1
			genes = deepcopy(chromosome1.genes)
			# Choose a random crossover position
			pos = np.random.random_integers(low=1, high=chromosome1.size - 2)
			# Crossover
			for i in range(pos):
				chromosome1.genes[i] = chromosome2.genes[i]
				chromosome2.genes[i] = genes[i]

	def mutate_chromosome(self, chromosome, prob_mutation, bounds):
		for i in range(chromosome.size):
			if (np.random.rand() <= prob_mutation):
				chromosome.genes[i] = chromosome.random_gene(bounds[i][0], bounds[i][1])
		
	def produce_offspring(self, bounds, prob_crossover, prob_mutation):
		# Check if the number of chromosomes is even
		even = True
		pairs = 0
		if (self.size % 2 == 0):
			pairs = self.size / 2
		else:
			pairs = (self.size + 1) / 2
			even = False
		# Select the pairs of parents and produce an offspring
		offspring = []
		for i in range(pairs):
			# Select parents via tournament selection
			idx1 = self.tournament_selection()
			idx2 = self.tournament_selection()
			chromosome1 = deepcopy(self.chromosomes[idx1])
			chromosome2 = deepcopy(self.chromosomes[idx2])
			# Crossover parents
			self.crossover_chromosomes(chromosome1, chromosome2, prob_crossover)
			# Mutate parents
			self.mutate_chromosome(chromosome1, prob_mutation, bounds)
			self.mutate_chromosome(chromosome2, prob_mutation, bounds)
			# Save new chromosomes
			offspring.append(chromosome1)
			offspring.append(chromosome2)
		if even:
			self.chromosomes = offspring
		else:
			self.chromosomes = offspring[:self.size-1]
			
	def score_chromosomes(self, indices, fixed, sim, exp, spinA, spinB, calc):
		func = partial(score_function, indices=indices, fixed=fixed, sim=sim, exp=exp, spinA=spinA, spinB=spinB, calc=calc)
		var = []
		for i in range(self.size): 
			var.append(self.chromosomes[i].genes)
		pool = Pool()
		score = pool.map(func, var)
		pool.close()
		pool.join()
		for i in range(self.size):
			self.chromosomes[i].score = score[i]	
				
	def sort_chromosomes(self):
		self.chromosomes.sort()