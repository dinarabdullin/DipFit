'''
Load the results of previous fitting
'''

import numpy as np
from supplement.constants import const


def chunkstring(string, length):
    return (string[0+i:length+i] for i in range(0, len(string), length))


def load_solution(filepath):
    solution = {}
    count = 0
    file = open(filepath + '/parameters.dat', 'r')
    # Skip a header
    next(file)
    for line in file:
        str = list(chunkstring(line, 16))
        parameter = {}
        name = const['variableNames'][count] 
        parameter['longname'] = str[0].strip()
        parameter['value'] = float(str[1]) * const['variableScales'][name]
        parameter['optimized'] = str[2].strip()
        parameter['precision'] = float(str[3]) * const['variableScales'][name]
        solution[name] = parameter
        count += 1
    return solution