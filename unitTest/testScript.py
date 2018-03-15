def arraytest(*npArray):
    import numpy as np
    print 'In testing script'
    shape = np.shape(npArray)
    print npArray
    print shape[1]
    
    a = np.arange(start = 0, stop = shape[1], step = 1.0, dtype='double')
    return a


def multiply(a,b):
    print "Will compute", a, "times", b
    return a*b
