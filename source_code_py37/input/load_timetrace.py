'''
Load an experimental RIDME time trace
'''

import numpy as np


def load_timetrace(filepath):
	xList = []
	yList = []
	file = open(filepath, 'r')
	for line in file:
		str = line.split()
		xList.append(float(str[0]))
		yList.append(float(str[1]))
	xExp = np.array(xList)	
	yExp = np.array(yList)	
	return [xExp, yExp]