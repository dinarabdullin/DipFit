'''
Genetic Algorithm: Saves the fit to the experimental spectrum
'''

def save_fit(xf, yf, xe, ye, filename):
    file = open(filename, 'w')
    if not (xe == []):
        file.write("{0:<12s} {1:<12s} {2:<12s} \n".format('f', 'exp', 'fit'))
        for i in range(xf.size):
            file.write('{0:<12.4f} {1:<12.4f} {2:<12.4f} \n'.format(xe[i], ye[i], yf[i]))
    else:
        file.write("{0:<12s} {1:<12s} \n".format('t', 'sim'))
        for i in range(xf.size):
            file.write('{0:<12.4f} {1:<12.4f} \n'.format(xf[i], yf[i]))
    file.close()