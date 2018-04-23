# This routine is required to return a vector of displacements from the orginal (t=0.0) location
# for each node in the list
# time - the non-dimensional time to set the movement for
# nodelist - list of indices for the nodes for which we need to move on the surface
# xyz - the unfiltered list of all nodes on this domain (index list x = xyz[inode*3 + 0]
# dxyz - the filtered list of all displacements needing to set (index list dx = dxyz[i*3 + 0]
def getBoundaryMovement(time, nodelist, xyz, dxyz):
    #Coords array(non-dimensional) [x,y,z]
    maxnodes = 0
    for i in range(0, maxnodes):
        inode = nodelist[i]
        x = xyz[inode*3 + 0]
        y = xyz[inode*3 + 1]
        z = xyz[inode*3 + 2]

        dxyz[i + 0] = 0.0
        dxyz[i + 1] = 0.0
        dxyz[i + 2] = 0.0
 
