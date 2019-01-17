'''
Genetic Algorithm: Saves the optimized values of fitting parameters
'''

def save_fitting_parameters(parameters, filename):    
    file = open(filename, 'w')
    file.write("{0:<16s} {1:<16s} {2:<16s} {3:<16s}\n".format('Parameter', 'Value', 'Optimized', 'Precision (+/-)'))
    for entry in parameters:
        file.write("{0:<16s} {1:<16.2f} {2:<16s} {3:<16.2f}\n".format(entry[0], entry[1], entry[2], entry[3]))
    file.close()