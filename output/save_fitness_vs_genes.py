'''
Genetic Algorithm: Saves the fitness as a function of individual genes
'''
	
from constants import const 
	
def save_fitness_vs_genes(fitness_vs_genes, gene_sets, directory):
	M = len(gene_sets)
	for i in range(M):
		# Number of genes in one set
		dim = len(gene_sets[i])
		# Create a new file
		filename = directory + 'parameter_errors'
		for j in range(dim):
			filename += '-' + gene_sets[i][j]
		filename += '.dat'
		file = open(filename, 'w')
		# Column names
		for j in range(dim):
			file.write("{0:<12s}".format(gene_sets[i][j]))
		file.write("{0:<12s}".format('RMSD'))
		file.write("\n")
		# Column content
		N = len(fitness_vs_genes[i]['fitness'])
		for k in range(N):
			for j in range(dim):
				name = gene_sets[i][j]
				par = fitness_vs_genes[i][name][k] / const['variableScales'][name]
				file.write("{0:<12.6f}".format(par))
			file.write("{0:<12.6f}".format(fitness_vs_genes[i]['fitness'][k]))
			file.write("\n")
		file.close()