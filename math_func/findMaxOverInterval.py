'''
Calculate root-mean-square-deviation (RMSD) between two signals
'''

import numpy as np

def findMaxOverInterval(y, x, xMin, xMax):
    ys = []
    for i in range(x.size):
        if (x[i] >= xMin) and (x[i] <= xMax):
            ys.append(y[i])
    ysMax = max(ys)	
    return ysMax