# This routine is required to return a vector of displacements from the orginal (t=0.0) location
# for each node in the list
# time - the non-dimensional time to set the movement for
# nodelist - list of indices for the nodes for which we need to move on the surface (indexes into xyz_base and dxyz)
# xyz_base - the unfiltered list of all nodes on this domain (index list x = xyz[inode*3 + 0] at time = 0
# xyz_current - the unfiltered list of all nodes on this domain (index list x = xyz[inode*3 + 0] at current time
# dxyz - the filtered list of all displacements needing to set (index list dx = dxyz[inode*3 + 0]
def getBoundaryMovement(time, nodelist, xyz_base, xyz_current, dxyz):
    #Coords array(non-dimensional) [x,y,z]
    maxnodes = 0
    for i in range(0, maxnodes):
        inode = nodelist[i]
        x = xyz_base[inode*3 + 0]
        y = xyz_base[inode*3 + 1]
        z = xyz_base[inode*3 + 2]

        dxyz[inode*3 + 0] = 0.0
        dxyz[inode*3 + 1] = 0.0
        dxyz[inode*3 + 2] = 0.0
 
