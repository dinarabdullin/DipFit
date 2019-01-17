'''
Calculate root-mean-square-deviation (RMSD) between two signals
'''

import numpy as np

def rmsd(sSim, sExp, f, fMin, fMax):
    rmsd = 0.0
    count = 0
    for i in range(sExp.size):
        if (np.absolute(f[i]) >= fMin) and (np.absolute(f[i]) <= fMax):
            rmsd += (sSim[i] - sExp[i])**2
            count += 1
    rmsd = np.sqrt(rmsd/float(count))	
    return rmsd