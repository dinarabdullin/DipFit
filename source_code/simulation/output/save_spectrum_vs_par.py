'''
Save the simulated spectra independence of some parameter
'''


def save_spectrum_vs_par(xs, par, ys, filename):
    file = open(filename, 'w')
    Nf = xs.size
    Npar = par.size
    file.write('{0:<12s}'.format(' '))
    for k in range(Npar):
        file.write('{0:<12.4f}'.format(par[k]))
    file.write('\n')
    for i in range(Nf):
        file.write('{0:<12.4f}'.format(xs[i]))
        for k in range(Npar):
            file.write('{0:<12.4f}'.format(ys[k][i]))
        file.write('\n')
    file.close()