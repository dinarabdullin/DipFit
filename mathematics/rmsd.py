'''
Calculate root-mean-square-deviation (RMSD) between two signals
'''

import numpy as np

def rmsd(ys, ye, x, xMin, xMax):
    rmsd = 0.0
    count = 0
    for i in range(ye.size):
        if (np.absolute(x[i]) >= xMin) and (np.absolute(x[i]) <= xMax):
            rmsd += (ys[i] - ye[i])**2
            count += 1
    rmsd = np.sqrt(rmsd/float(count))	
    return rmsd