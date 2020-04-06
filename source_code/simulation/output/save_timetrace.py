'''
Save the experimental and simulated time traces
'''


def save_timetrace(t_sim, sig_sim, t_exp, sig_exp, filename):
    file = open(filename, 'w')
    if not (t_exp == []):
        file.write("{0:<12s} {1:<12s} {2:<12s} \n".format('t', 'exp', 'sim'))
        for i in range(t_sim.size):
            file.write('{0:<12.4f} {1:<12.4f} {2:<12.4f} \n'.format(t_exp[i], sig_exp[i], sig_sim[i]))
    else:
        file.write("{0:<12s} {1:<12s} \n".format('t', 'sim'))
        for i in range(t_sim.size):
            file.write('{0:<12.4f} {1:<12.4f} \n'.format(t_sim[i], sig_sim[i]))
    file.close()