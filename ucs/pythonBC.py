def pythonBC(QL, QR, neqn, naux, CoordsXYZ):
    import numpy as np

    # You have two options with setting BCs in this way
    # 1) You can set QL which is physically ON the wall at CoordsXYZ - hard set BC
    #  --- This is used often for isothermal heat transfer BCs
    # 2) You can set QR which is a ghost cell to control the wall flux - soft set BC
    #  --- This is used often for farfield characteristic type BCs

    print neqn + naux
    qinfSize = neqn + naux
    print 'Vector size', qinfSize

    #Coords array(non-dimensional) [x,y,z]
    x = CoordsXYZ[0]
    y = CoordsXYZ[1]
    z = CoordsXYZ[2]
    print x,y,z

    #Note that computeAuxiliaryVariables() is called after the set to ensure consistency
    # for all auxiliary values, we only have to set the conservative variables appropriately
    
    QR[0] = 1.5
