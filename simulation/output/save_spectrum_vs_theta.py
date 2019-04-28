'''
Saves the dependence of the simulated spectrum on the angle 'theta'
'''

def save_spectrum_vs_theta(xs, theta, ys_theta, filename):
    file = open(filename, 'w')
    Nf = xs.size
    Ntheta = theta.size
    file.write('{0:<12s}'.format(' '))
    for k in range(Ntheta):
        file.write('{0:<12.4f}'.format(theta[k]))
    file.write('\n')
    for i in range(Nf):
        file.write('{0:<12.4f}'.format(xs[i]))
        for k in range(Ntheta):
            file.write('{0:<12.4f}'.format(ys_theta[k][i]))
        file.write('\n')
    file.close()