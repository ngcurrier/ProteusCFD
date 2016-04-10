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
    print 'Reading ' + str(int(nhex)) + ' hexes'

    triCount = 0
    quadCount = 0
    tetCount = 0
    pyramidCount = 0
    prismCount = 0
    hexCount = 0

    elements = []
    
    i = 0;
    while i < len(elemData):
        if elemData[i] == TRI:
            i = i + 3 + 1
            elements.append(Tri())
            elements[-1].nodes = elemData[i+1:i+4]
            triCount = triCount + 1
            continue
        elif elemData[i] == QUAD:
            i = i + 4 +1
            elements.append(Quad())
            elements[-1].nodes = elemData[i+1:i+5]
            quadCount = quadCount + 1
            continue
        elif elemData[i] == TET:
            i = i + 4 +1
            elements.append(Tet())
            elements[-1].nodes = elemData[i+1:i+5]
            tetCount = tetCount + 1
            continue
        elif elemData[i] == PYRAMID:
            i = i + 5 + 1
            elements.append(Pyramid())
            elements[-1].nodes = elemData[i+1:i+6]
            pyramidCount = pyramidCount + 1
            continue
        elif elemData[i] == PRISM:
            i = i + 6 + 1
            elements.append(Prism())
            elements[-1].nodes = elemData[i+1:i+7]
            prismCount = prismCount + 1
            continue
        elif elemData[i] == HEX:
            i = i + 8 + 1
            elements.append(Hex())
            elements[-1].nodes = elemData[i+1:i+9]
            hexCount = hexCount + 1
            continue

    print 'Found ' + str(triCount) + ' Tri'
    print 'Found ' + str(quadCount) + ' Quad'
    print 'Found: ' + str(tetCount) + ' Tet'
    print 'Found: ' + str(pyramidCount) + ' Pyramid'
    print 'Found: ' + str(prismCount) + ' Prism'
    print 'Found: ' + str(hexCount) + ' Hex'

    #for elm in elements:
    #    print elm.getName()
    #    print elm.nodes
    
    
class Elem():
    def __init__(self):
        self.nnodes = 0
        self.nodes = []

    def getNnodes(self):
        return self.nnodes;

    def getName(self):
        raise NotImplementedError

class Tri(Elem):
    def __init__(self):
        self.nnodes = 3
        self.nodes = [0,0,0]

    def getName(self):
        return 'tri'

class Quad(Elem):
    def __init__(self):
        self.nnodes = 4;
        self.nodes = [0,0,0,0]

    def getName(self):
        return 'quad'
        
class Tet(Elem):
    def __init__(self):
        self.nnodes = 4
        self.nodes = [0,0,0,0]

class Pyramid(Elem):
    def __init__(self):
        self.nnodes = 5
        self.nodes = [0,0,0,0,0]

class Prism(Elem):
    def __init__(self):
        self.nnodes = 6
        self.nodes  = [0,0,0,0,0,0]

class Hex(Elem):
    def __init__(self):
        self.nnodes = 8
        self.nodes = [0,0,0,0,0,0,0,0]

    def getName(self):
        return 'hex'
        
def main():
    loadHDF5File("bump.0.h5")
    
if __name__ == "__main__":
    main()
