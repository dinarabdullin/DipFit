'''
Genetic Algorithm: Saves the fitness as a function of the number of generations
'''
	
def save_fitness(fitness, filename):
    file = open(filename, 'w')
    for i in range(fitness.size):
        file.write('{0:<12d} {1:<12.6f} \n'.format(i+1, fitness[i]))
    file.close()