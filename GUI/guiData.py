import h5reader
import vtk
from defines import *

class GuiData():
    def __init__(self):
        self.params = []
        self.bcTypes = {"velocity", "pressure", "farfield", "symmetry"}
        self.elemGroups = []
        self.mesh = None
        self.meshes = []

    def loadMesh(self, path, casename):
        self.np = h5reader.checkNumberPartitions(path, casename)
        print 'Number of processors: ' + str(self.np)
        for i in range(0, self.np):
            self.meshes.append(h5reader.loadHDF5File(path, casename, i))
        self.mesh  = self.concatenateMeshes()
        # TODO: free memory here once concatenation is done
        self.elemGroups = self.mesh.elementGroups 

    def concatenateMeshes(self):
        m = h5reader.Mesh()
        nodeNumOffset = 0 
        m.nelem = self.meshes[0].nelem # global counter
        for i in range(0, self.np):
            m.nnodes += self.meshes[i].nnodes
            m.ntri += self.meshes[i].ntri
            m.nquad += self.meshes[i].nquad
            m.ntet += self.meshes[i].ntet
            m.npyramid += self.meshes[i].npyramid
            m.nprism += self.meshes[i].nprism
            m.nhex += self.meshes[i].nhex
            print 'Concatenating mesh: ' + str(i)
            for elem in self.meshes[i].elements:
                etemp = elem
                # shift node numbering so we still reference the same node in the now global list
                nodesShifted = []
                for node in etemp.nodes:
                    nodesShifted.append(int(node)+nodeNumOffset)
                if len(nodesShifted) == 0: continue
                etemp.setNodes(nodesShifted)    
                m.elements.append(etemp)
            for coord in self.meshes[i].coords:
                m.coords.append(coord)
            print 'Size of coords list: ' + str(len(m.coords))
            for tag in self.meshes[i].factags:
                m.factags.append(tag)
            nodeNumOffset += self.meshes[i].nnodes
        m.elementGroups = h5reader.uniqify(m.factags)
        m.elementGroups = sorted(m.elementGroups)
        print 'Unique element groups -- (+) boundaries, (-) volumes: ' + str(m.elementGroups)
        # count all positive factags
        m.nfactags = 0
        for item in m.elementGroups:
            if int(item) > 0:
                m.nfactags = m.nfactags + 1

        # check that total global elements matches with what we have
        totalElem = (m.ntri + m.nquad + m.ntet + m.npyramid + m.nprism + m.nhex)
        print 'Global number of elements expected ' + str(m.nelem) + ' number found is ' + str(totalElem)
        print 'Global factags read ' + str(len(m.factags))
        print 'Total duplicates present are: ' + str(totalElem - m.nelem)

        print 'Check coordinate number: ' + str(len(m.coords)/3) + ' nnodes:' + str(m.nnodes)
        
        
        # Get max x value
        maxx = m.coords[0]
        minx = m.coords[0]
        for i in range(0,m.nnodes):
            maxx = max(maxx, m.coords[i*3 + 0])
            minx = min(minx, m.coords[i*3 + 0])
        # Get max y value
        maxy = m.coords[1]
        miny = m.coords[1]
        for i in range(0,m.nnodes):
            maxy = max(maxy, m.coords[i*3 + 1])
            miny = min(miny, m.coords[i*3 + 1])
        # Get max z value
        maxz = m.coords[2]
        minz = m.coords[2]
        for i in range(0,m.nnodes):
            maxz = max(maxz, m.coords[i*3 + 2])
            minx = min(minz, m.coords[i*3 + 2])
                       
        print 'Xrange: ' + str(minx) + ':' + str(maxx)
        print 'Yrange: ' + str(miny) + ':' + str(maxy)
        print 'Zrange: ' + str(minz) + ':' + str(maxz)
                                              
        return m
            
    #remove data contained re: mesh
    def clear(self):
        self.mesh = None
        self.meshes = []
        self.elemGroups = []

    #builds a set of vtk grid objects based on factag defined
    #element groups, i.e. boundaries and volumes
    def buildVTKGrids(self):
        grids = []
        for tag in self.mesh.elementGroups:
            grids.append(self.buildPartialVTKGrid(tag))
        #grids.append(self.buildFullVTKGrid())
        return grids
        
    def getPointsWithFactag(self, factag):
        #loop over elements, add unique points to list
        seen = {}
        result = []
        for element in self.getElementsWithFactag(factag):
            for node in element.nodes:
                if int(node) in seen: continue
                seen[int(node)] = 1
                result.append(int(node))
        return result

    def getElementsWithFactag(self, factag):
        #loop over elements, add unique elements to list
        result = []
        for element in self.mesh.elements:
            if element.factag == factag:
                result.append(element)
        return result


    # Builds a VTK grid based on factags
    def buildPartialVTKGrid(self, factag):
        #get points required for factag
        pointsList = self.getPointsWithFactag(factag)

        #create a lookup table so we can map the
        #cells from the global list to a local list
        points = vtk.vtkPoints()
        localIdx = 0
        ptMap = {}
        for pt in pointsList:
            ptMap[int(pt)] = localIdx
            localIdx = localIdx + 1
            p = self.mesh.coords[(pt*3):(pt*3+3)]
            points.InsertNextPoint(p)
        
        vtkgrid = vtk.vtkUnstructuredGrid()
        vtkgrid.SetPoints(points)

        #get elements that have desired factag
        felements = self.getElementsWithFactag(factag)

        #build the vtk elements
        for element in felements:
            type = element.getType()
            nodes = element.nodes
            if type == eTypes.TRI:
                cell = vtk.vtkTriangle()
            elif type == eTypes.QUAD:
                cell = vtk.vtkQuad()
            elif type == eTypes.TET:
                cell = vtk.vtkTetra()
            elif type == eTypes.PYRAMID:
                cell = vtk.vtkPyramid()
            elif type == eTypes.PRISM:
                cell = vtk.vtkWedge()  #prism
            elif type == eTypes.HEX:
                cell = vtk.vtkHexahedron()
            else:
                raise # throw an exception
            j = 0
            for n in nodes:
                localId = ptMap[int(n)]
                cell.GetPointIds().SetId(j,localId)
                j = j+1
            vtkgrid.InsertNextCell(cell.GetCellType(), cell.GetPointIds())
        return vtkgrid

    # Builds a VTK grid as a whole, no factag piecewise build
    def buildFullVTKGrid(self):
        # Create the points for VTK
        points = vtk.vtkPoints()
        for i in range(0, len(self.mesh.coords)/3):
            p = self.mesh.coords[(i*3):(i*3+3)]
            points.InsertNextPoint(p)

        #add the points and cells to unstructured grid
        vtkgrid = vtk.vtkUnstructuredGrid()
        vtkgrid.SetPoints(points)

        #add the VTK elements to the mesh
        for element in self.mesh.elements:
            type = element.getType()
            nodes = element.nodes
            cell = vtk.vtkTriangle()
            if type == eTypes.TRI:
                cell = vtk.vtkTriangle()
            elif type == eTypes.QUAD:
                cell = vtk.vtkQuad()
            elif type == eTypes.TET:
                cell = vtk.vtkTetra()
            elif type == eTypes.PYRAMID:
                cell = vtk.vtkPyramid()
            elif type == eTypes.PRISM:
                cell = vtk.vtkWedge()  #prism
            elif type == eTypes.HEX:
                cell = vtk.vtkHexahedron()
            else:
                raise # throw an exception
            j = 0
            for n in nodes:
                cell.GetPointIds().SetId(j,n)
                j = j+1
            vtkgrid.InsertNextCell(cell.GetCellType(), cell.GetPointIds())
        return vtkgrid
