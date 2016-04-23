#!/usr/bin/python

import h5py as h5
import numpy as np
import vtk
from defines import *

#todo: glob all files with correct filename pattern for parallel read

def uniqify(seq): 
   seen = {}
   result = []
   for i in range(0,len(seq)):
       marker = int(seq[i])
       # in old Python versions:
       # if seen.has_key(marker)
       # but in new ones:
       if marker in seen: continue
       seen[marker] = 1
       result.append(int(seq[i]))
   return result
        
def loadHDF5File(casename):
    filename = casename + '.0.h5'
    h5f = h5.File(filename, "r")

    m = Mesh();

    m.nelem = h5f['Global Number Of Elements'][:]
    elemData = h5f['Mesh/Element Data'][:]
    m.nfactags = h5f['Number Of Factags'][:]
    m.factags = h5f['Mesh/Element Factags'][:]
    m.coords = h5f['Mesh/Nodal Coordinates'][:]

    m.nprism = h5f['Mesh/Number Of Prisms'][:]
    m.npyramid = h5f['Mesh/Number Of Pyramids'][:]
    m.nquad = h5f['Mesh/Number Of Quadrilaterals'][:]
    m.ntet = h5f['Mesh/Number Of Tetrahedron'][:]
    m.ntri = h5f['Mesh/Number Of Prisms'][:]
    m.nhex = h5f['Mesh/Number Of Hexahedron'][:]
    
    m.printElemCounts()
    
    triCount = 0
    quadCount = 0
    tetCount = 0
    pyramidCount = 0
    prismCount = 0
    hexCount = 0

    elements = []

    i = 0;
    icell = 0;
    while i < len(elemData):
        if elemData[i] == eTypes.TRI:
            i = i + 3 + 1
            m.elements.append(Tri())
            m.elements[-1].nodes = elemData[i+1:i+4]
            m.elements[-1].factag = m.factags[icell]
            triCount = triCount + 1
            icell = icell + 1
            continue
        elif elemData[i] == eTypes.QUAD:
            i = i + 4 +1
            m.elements.append(Quad())
            m.elements[-1].nodes = elemData[i+1:i+5]
            m.elements[-1].factag = m.factags[icell]
            quadCount = quadCount + 1
            icell = icell + 1
            continue
        elif elemData[i] == eTypes.TET:
            i = i + 4 +1
            m.elements.append(Tet())
            m.elements[-1].nodes = elemData[i+1:i+5]
            m.elements[-1].factag = m.factags[icell]
            tetCount = tetCount + 1
            icell = icell + 1
            continue
        elif elemData[i] == eTypes.PYRAMID:
            i = i + 5 + 1
            m.elements.append(Pyramid())
            m.elements[-1].nodes = elemData[i+1:i+6]
            m.elements[-1].factag = m.factags[icell]
            pyramidCount = pyramidCount + 1
            icell = icell + 1
            continue
        elif elemData[i] == eTypes.PRISM:
            i = i + 6 + 1
            m.elements.append(Prism())
            m.elements[-1].nodes = elemData[i+1:i+7]
            m.elements[-1].factag = m.factags[icell]
            prismCount = prismCount + 1
            icell = icell + 1
            continue
        elif elemData[i] == eTypes.HEX:
            i = i + 8 + 1
            m.elements.append(Hex())
            m.elements[-1].nodes = elemData[i+1:i+9]
            m.elements[-1].factag = m.factags[icell]
            hexCount = hexCount + 1
            icell = icell + 1
            continue
        else:
            raise

    print 'Unique element groups -- (+) boundaries, (-) volumes: ' + str(m.nfactags)
    m.elementGroups = uniqify(m.factags);
    print m.elementGroups
    
    print 'Found ' + str(triCount) + ' Tri'
    print 'Found ' + str(quadCount) + ' Quad'
    print 'Found: ' + str(tetCount) + ' Tet'
    print 'Found: ' + str(pyramidCount) + ' Pyramid'
    print 'Found: ' + str(prismCount) + ' Prism'
    print 'Found: ' + str(hexCount) + ' Hex'

    #for elm in elements:
    #    print elm.getName()
    #    print elm.nodes
    return m

class Mesh():
    
    def __init__(self):
        self.nnodes = 0
        self.nfactags = 0
        self.nelem = 0
        self.ntri = 0
        self.nquad = 0
        self.ntet = 0
        self.npyramid = 0
        self.nprism = 0
        self.nhex = 0
        self.elements = []
        self.coords = []
        self.factags = []
        self.elementGroups = []
        
    def printElemCounts(self):
        print 'Reading ' + str(int(self.nprism)) + ' prisms'
        print 'Reading ' + str(int(self.npyramid)) + ' pyramids'
        print 'Reading ' + str(int(self.nquad)) + ' quads'
        print 'Reading ' + str(int(self.ntet)) + ' tets'
        print 'Reading ' + str(int(self.ntri)) + ' triangles'
        print 'Reading ' + str(int(self.nhex)) + ' hexes'

            
class Elem():
    def __init__(self):
        self.nnodes = 0
        self.factag = -999
        self.nodes = []

    def getNnodes(self):
        return self.nnodes;

    def getFactag(self):
        return factag

    def getName(self):
        raise NotImplementedError

    def getType(self):
        raise NotImplementedError

class Tri(Elem):
    def __init__(self):
        self.nnodes = 3
        self.nodes = [0,0,0]

    def getName(self):
        return 'tri'

    def getType(self):
        return eTypes.TRI
 
class Quad(Elem):
    def __init__(self):
        self.nnodes = 4;
        self.nodes = [0,0,0,0]

    def getName(self):
        return 'quad'
        
    def getType(self):
        return eTypes.QUAD

class Tet(Elem):
    def __init__(self):
        self.nnodes = 4
        self.nodes = [0,0,0,0]

    def getName(self):
        return 'tet'
        
    def getType(self):
        return eTypes.TET

class Pyramid(Elem):
    def __init__(self):
        self.nnodes = 5
        self.nodes = [0,0,0,0,0]

    def getName(self):
        return 'pyramid'
        
    def getType(self):
        return eTypes.PYRAMID

class Prism(Elem):
    def __init__(self):
        self.nnodes = 6
        self.nodes  = [0,0,0,0,0,0]

    def getName(self):
        return 'prism'
    
    def getType(self):
        return eTypes.PRISM


class Hex(Elem):
    def __init__(self):
        self.nnodes = 8
        self.nodes = [0,0,0,0,0,0,0,0]

    def getName(self):
        return 'hex'
        
    def getType(self):
        return eTypes.HEX

def main():
    loadHDF5File("bump")
    
if __name__ == "__main__":
    main()
