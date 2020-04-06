'''
Save the simulated and experimental spectra
'''


def save_spectrum(f_sim, spc_sim, f_exp, spc_exp, filename):
    file = open(filename, 'w')
    if not (f_exp == []):
        file.write("{0:<12s} {1:<12s} {2:<12s} \n".format('f', 'exp', 'sim'))
        for i in range(f_sim.size):
            file.write('{0:<12.4f} {1:<12.4f} {2:<12.4f} \n'.format(f_exp[i], spc_exp[i], spc_sim[i]))
    else:
        file.write("{0:<12s} {1:<12s} \n".format('t', 'sim'))
        for i in range(f_sim.size):
            file.write('{0:<12.4f} {1:<12.4f} \n'.format(f_sim[i], spc_sim[i]))
    file.close()