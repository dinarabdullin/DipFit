'''
Save the score in dependence of individual parameters
'''
	
from supplement.constants import const 

	
def save_score_vs_par(score_vs_par, par, directory):
	M = len(par)
	for i in range(M):
		# Number of genes in one set
		dim = len(par[i])
		# Create a new file
		filename = directory + 'parameter_errors'
		for j in range(dim):
			filename += '-' + par[i][j]
		filename += '.dat'
		file = open(filename, 'w')
		# Column names
		for j in range(dim):
			file.write("{0:<12s}".format(par[i][j]))
		file.write("{0:<12s}".format('RMSD'))
		file.write("\n")
		# Column content
		N = len(score_vs_par[i]['score'])
		for k in range(N):
			for j in range(dim):
				name = par[i][j]
				val = score_vs_par[i][name][k] / const['variableScales'][name]
				file.write("{0:<12.6f}".format(val))
			file.write("{0:<12.6f}".format(score_vs_par[i]['score'][k]))
			file.write("\n")
		file.close()