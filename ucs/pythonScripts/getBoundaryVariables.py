# Set the boundary conditions at the surface, user may choose to set QL, QR, or both
# QL - raw vector on the surface (hard set BC)
# QR - raw vector on the ghost edges (soft set BC)
# neqn - number of conservative equations
# naux - number of auxiliary (helper) variables
# wallXYZ - X,Y,Z location of the wall node
# wallAvec - normal vector of face at the wall (nx, ny, nz, area)
def getBoundaryVariables(QL, QR, neqn, naux, wallXYZ, wallAvec):
    import numpy as np

    print neqn+ naux
    qinfSize = neqn + naux
    print 'Vector size', qinfSize

    # QL is the raw value the raw boundary variables (non-dimensional) - hard set BC
    # QR is the value in the ghost cell - soft set BC
    # computeAuxiliaryVariables() is called after return, don't bother setting aux vars
    QR[0] = 0.00
    QL[0] = 0.00
    
    #Coords array(non-dimensional) [x,y,z]
    x = wallXYZ[0]
    y = wallXYZ[1]
    z = wallXYZ[2]
    print x,y,z

    # wall avec is the area vector (normal & area for the face)
    nx = wallAvec[0]
    ny = wallAvec[1]
    nz = wallAvec[2]
    a = wallAvec[3]

