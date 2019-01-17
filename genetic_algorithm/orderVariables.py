'''
Make a list with all variables in a pre-defined order
'''

import numpy as np
from constants import const

def order_variables(var, indices, fixed):
    values = []
    optimized = []
	for name in const['variableNames']:
		index = indices[name]
		if not (index == -1):
			values.append(var[index] / const['variableScales'])
			optimized.append('Y')
		else:
			values.append(fixed[name] / const['variableScales'])
			optimized.append('N')
    var_ordered = []
    for i in range(len(indices)):
        var_ordered.append([const['longVariableNames'][i], values[i], optimized[i]])
    return var_ordered