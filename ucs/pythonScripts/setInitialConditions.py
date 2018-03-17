def setInitialConditions(Qinf, Neqn, Qset, CoordsXYZ):
    import numpy as np
    
    #Qinf is the raw value of farfield variables (non-dimensional) - DO NOT MODIFY
    qinfSize = np.shape(Qinf)[1]
    
    #Coords array(non-dimensional) [x,y,z]
    x = Coords[0]
    y = Coords[1]
    z = Coords[2]

    #Note that computeAuxiliaryVariables() is called after the set to ensure consistency
    #You do not have to set those locations to anything meaningfu
    naux = qinfSize - Neqn

    #Qset should be set to desired values (non-dimensional)
    Qset[0] = 10.0
