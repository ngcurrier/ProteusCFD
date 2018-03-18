def setInitialConditions(Qinf, neqn, naux, Qset, CoordsXYZ):
    import numpy as np

    print Qinf[0]
    
    #Qinf is the raw value of farfield variables (non-dimensional) - DO NOT MODIFY
    print neqn+ naux
    qinfSize = neqn + naux
    print 'Vector size', qinfSize
    
    #Coords array(non-dimensional) [x,y,z]
    x = CoordsXYZ[0]
    y = CoordsXYZ[1]
    z = CoordsXYZ[2]
    print x,y,z

    #Note that computeAuxiliaryVariables() is called after the set to ensure consistency
    #You do not have to set those locations to anything meaningfu
    naux = qinfSize - neqn

    #Qset should be set to desired values (non-dimensional)
    Qset[0] = 10.0

    return 1
