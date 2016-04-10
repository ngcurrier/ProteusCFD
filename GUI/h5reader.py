#!/usr/bin/python

import h5py as h5
import numpy as np

#define element types
TRI = 0
QUAD = 1
TET = 2
PYRAMID = 3
PRISM = 4
HEX = 5

def loadHDF5File(filename):
    h5f = h5.File(filename, "r")
    nelem = h5f['Global Number Of Elements'][:]
    elemData = h5f['Mesh/Element Data'][:]
    factags = h5f['Mesh/Element Factags'][:]
    coords = h5f['Mesh/Nodal Coordinates'][:]

    nprism = h5f['Mesh/Number Of Prisms'][:]
    npyramid = h5f['Mesh/Number Of Pyramids'][:]
    nquad = h5f['Mesh/Number Of Quadrilaterals'][:]
    ntet = h5f['Mesh/Number Of Tetrahedron'][:]
    ntri = h5f['Mesh/Number Of Prisms'][:]
    nhex = h5f['Mesh/Number Of Hexahedron'][:]
    
    print 'Reading ' + str(int(nprism)) + ' prisms'
    print 'Reading ' + str(int(npyramid)) + ' pyramids'
    print 'Reading ' + str(int(nquad)) + ' quads'
    print 'Reading ' + str(int(ntet)) + ' tets'
    print 'Reading ' + str(int(ntri)) + ' triangles'
    print 'Reading ' + str(int(nhex)) + ' nhexes'

    elements = [tri(), quad(), hex()]

    for elm in elements:
        print elm.getName()
        print elm.nodes
    
    tri1 = tri()
    print 'tri nodes: ' + str(tri1.getNnodes())
    
class elem():
    def __init__(self):
        self.nnodes = 0
        self.nodes = []

    def getNnodes(self):
        return self.nnodes;

    def getName(self):
        raise NotImplementedError

class tri(elem):
    def __init__(self):
        self.nnodes = 3
        self.nodes = [0,0,0]

    def getName(self):
        return 'tri'
        

class quad(elem):
    def __init__(self):
        self.nnodes = 4;
        self.nodes = [0,0,0,0]

    def getName(self):
        return 'quad'
        
class tet(elem):
    def __init__(self):
        self.nnodes = 4
        self.nodes = [0,0,0,0]

class pyramid(elem):
    def __init__(self):
        self.nnodes = 5
        self.nodes = [0,0,0,0,0]

class prism(elem):
    def __init__(self):
        self.nnodes = 6
        self.nodes  = [0,0,0,0,0,0]

class hex(elem):
    def __init__(self):
        self.nnodes = 8
        self.nodes = [0,0,0,0,0,0,0,0]

    def getName(self):
        return 'hex'
        
def main():
    loadHDF5File("bump.0.h5")
    
if __name__ == "__main__":
    main()
