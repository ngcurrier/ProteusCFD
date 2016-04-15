import h5reader
import vtk
from defines import *

class GuiData():
    def __init__(self):
        self.params = []
        self.bcTypes = []

    def loadMesh(self):
        self.caseName = 'bump'
        self.mesh = h5reader.loadHDF5File(self.caseName)
        self.buildVTKGrid();
        
    def buildVTKGrid(self):
        # Create the points for VTK
        points = vtk.vtkPoints()
        for i in range(0, len(self.mesh.coords)/3):
            p = self.mesh.coords[(i*3):(i*3+3)]
            points.InsertNextPoint(p)

        #add the points and cells to unstructured grid
        self.vtkgrid = vtk.vtkUnstructuredGrid()
        self.vtkgrid.SetPoints(points)

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
            self.vtkgrid.InsertNextCell(cell.GetCellType(), cell.GetPointIds())

    def getVTKGrid(self):
        return self.vtkgrid
