'''
Genetic Algorithm: Saves the score as a function of the number of generations
'''
	
def save_score(score, filename):
	file = open(filename, 'w')
	for i in range(score.size):
		file.write('{0:<12d} {1:<12.6f} \n'.format(i+1, score[i]))
	file.close()