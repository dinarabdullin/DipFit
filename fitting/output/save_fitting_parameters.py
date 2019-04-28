'''
Genetic Algorithm: Saves the optimized values of fitting parameters
'''

from supplement.constants import const

def save_fitting_parameters(parameters, filename):    
	file = open(filename, 'w')
	file.write("{0:<16s} {1:<16s} {2:<16s} {3:<16s}\n".format('Parameter', 'Value', 'Optimized', 'Precision (+/-)'))
	for name in const['variableNames']:
		parameter = parameters[name]
		file.write("{0:<16s} {1:<16.3f} {2:<16s} {3:<16.3f}\n".format(parameter['longname'], parameter['value']/const['variableScales'][name], parameter['optimized'], parameter['precision']/const['variableScales'][name]))
	file.close()