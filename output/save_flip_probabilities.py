'''
Saves the flip probability vs g-factor and temperature
'''
    
def save_flip_probabilities(g, l, T, filename):
	file = open(filename, 'w')
	Nt = len(T)
	Ng = g.size
	# Column names
	file.write("{0:<12s}".format('g'))
	for i in range(Nt):
		label = str(T[i]) + ' K'
		file.write("{0:<12s}".format(label))
	file.write("\n")
	# Column content
	for j in range(Ng):
		file.write("{0:<12.4f}".format(g[j]))
		for i in range(Nt):
			file.write("{0:<12.4f}".format(l[i][j]))
		file.write("\n")
	file.close()