'''
Checks whether a function y(x) has symmetric boundaries
'''

def symmetric_boundaries(x, y):
    xsym = []
    ysim = []
    if (x[0] == -x[-1]):
        xsym = x
        ysym = y
    if (x[0] < -x[0]):
        xsym = x[1:]
        ysym = y[1:]
    if (x[0] > -x[0]):
        xsym = x[:-1]
        ysym = y[:-1]    
    return [xsym, ysym]