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

def setInitialConditions(Qinf, neqn, naux, Qset, CoordsXYZ):
    
    #Qinf is the raw value of farfield variables (non-dimensional) - DO NOT MODIFY
    qinfSize = np.shape(Qinf)[1]
    print 'Vector size' + str(qinfsize)
    
    #Coords array(non-dimensional) [x,y,z]
    x = Coords[0]
    y = Coords[1]
    z = Coords[2]
    print x,y,z

    #Note that computeAuxiliaryVariables() is called after the set to ensure consistency
    #You do not have to set those locations to anything meaningfu
    naux = qinfSize - neqn

    #Qset should be set to desired values (non-dimensional)
    Qset[0] = 10.0
