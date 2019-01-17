'''
Saves the simulated RIDME time trace
'''

def save_timetrace(xs, ys, xe, ye, filename):
    file = open(filename, 'w')
    if not (xe == []):
        file.write("{0:<12s} {1:<12s} {2:<12s} \n".format('t', 'exp', 'sim'))
        for i in range(xs.size):
            file.write('{0:<12.4f} {1:<12.4f} {2:<12.4f} \n'.format(xe[i], ye[i], ys[i]))
    else:
        file.write("{0:<12s} {1:<12s} \n".format('t', 'sim'))
        for i in range(xs.size):
            file.write('{0:<12.4f} {1:<12.4f} \n'.format(xs[i], ys[i]))
    file.close()