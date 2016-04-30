import h5reader
import vtk
from defines import *

class GuiData():
    def __init__(self):
        self.params = []
        self.bcTypes = {"velocity", "pressure", "farfield", "symmetry"}
        self.elemGroups = []

    def loadMesh(self, path, casename):
        self.mesh = h5reader.loadHDF5File(path, casename)
        self.elemGroups = self.mesh.elementGroups

    #remove data contained re: mesh
    def clear(self):
        self.mesh = None
        self.elemGroups = []

    #builds a set of vtk grid objects based on factag defined
    #element groups, i.e. boundaries and volumes
    def buildVTKGrids(self):
        grids = []
        for tag in self.mesh.elementGroups:
            grids.append(self.buildPartialVTKGrid(tag))
        return grids
        
    def getPointsWithFactag(self, factag):
        #loop over elements, add unique points to list
        seen = {}
        result = []
        for element in self.mesh.elements:
            nodes = element.nodes
            efactag = element.factag
            if factag != efactag: continue
            for node in nodes:
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


    def buildPartialVTKGrid(self, factag):
        #get points required for factag
        pointsList = self.getPointsWithFactag(factag)

        #create a lookup table so we can map the
        #cells from the global list to a local list
        localIdx = 0
        ptMap = {}
        for pt in pointsList:
            ptMap[int(pt)] = localIdx
            localIdx = localIdx + 1
        
        vtkgrid = vtk.vtkUnstructuredGrid()

        #create the points for VTK
        points = vtk.vtkPoints()
        for i in range(0, len(pointsList)):
            ptid = pointsList[i]
            p = self.mesh.coords[(ptid*3):(ptid*3+3)]
            points.InsertNextPoint(p)
        vtkgrid.SetPoints(points)

        #get elements that have desired factag
        felements = self.getElementsWithFactag(factag)

        #build the vtk elements
        for element in felements:
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
                localId = ptMap[int(n)]
                cell.GetPointIds().SetId(j,localId)
                j = j+1
            vtkgrid.InsertNextCell(cell.GetCellType(), cell.GetPointIds())
        return vtkgrid
        
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
