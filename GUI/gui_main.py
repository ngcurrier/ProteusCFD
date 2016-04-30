#!/usr/bin/env python
 
import sys
import pickle # we use pickle to save run information, etc.
import vtk
from guiData import *
from BCTab import BCTab
from PyQt4 import QtCore, QtGui
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from defines import *
import functools


class MainWindow(QtGui.QMainWindow):
 
    def __init__(self, parent = None):
        QtGui.QMainWindow.__init__(self, parent)

        self.state = {"test", "stuff", "other"}
        self.data = GuiData()
         
        self.frame = QtGui.QFrame()
        self.resize(1000,600)
        self.setWindowTitle('ProteusCFD')
        self.move(100,100)

        self.vl = QtGui.QGridLayout()
        self.vl.setSpacing(10)

        self.frame.setLayout(self.vl)
        self.setCentralWidget(self.frame)

        #create tabs dialog and add it to the layout
        self.tabs = QtGui.QTabWidget()
        self.vl.addWidget(self.tabs,0,0)

        #Setup our tree view
        self.treeWidget = QtGui.QTreeWidget()
        self.treeWidget.setHeaderHidden(True)
        self.addItems(self.treeWidget.invisibleRootItem())
        self.treeWidget.itemChanged.connect (self.handleChanged)
        self.tabs.addTab(self.treeWidget, "Physics Options")

        #setup solution control tab
        SolutionControlTab = QtGui.QWidget()
        self.tabs.addTab(SolutionControlTab, "Solution Control")

        #setup boundary condition tab
        BCTab1 = BCTab(self.data, self.visibilityChanged)
        self.tabs.addTab(BCTab1, "Boundary Conditions")
        
        #draw the tabs
        self.tabs.show()

        #Setup our menu and actions
        exitAction = QtGui.QAction('Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(self.close)

        saveAction = QtGui.QAction('Save', self)
        saveAction.setShortcut('Ctril+S')
        saveAction.setStatusTip('Save config')
        saveAction.triggered.connect(self.save)

        loadAction = QtGui.QAction('Load', self)
        loadAction.setStatusTip('Load config')
        loadAction.triggered.connect(self.load)

        openFileAction = QtGui.QAction('Import', self)
        openFileAction.setStatusTip('Import Grid')
        openFileAction.triggered.connect(self.selectFile)
        
        self.menubar = self.menuBar()
        fileMenu = self.menubar.addMenu('&File')
        fileMenu.addAction(exitAction)
        fileMenu.addAction(saveAction)
        fileMenu.addAction(loadAction)
        fileMenu.addAction(openFileAction)
        
        self.statusBar().showMessage('Waiting...')
        
        #Setup our vtk widget
        self.vtkWidget = QVTKRenderWindowInteractor(self.frame)
        self.vl.addWidget(self.vtkWidget, 0, 300)

        self.ren = vtk.vtkRenderer()
        self.ren.SetBackground(.5,.8,.5)       
        self.ren.ResetCamera()

        self.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
        self.iren = self.vtkWidget.GetRenderWindow().GetInteractor()

        #read the mesh
        meshFile = "bump.0.h5"
        parts = meshFile.split('.')
        self.casename = parts[0]
        self.data.loadMesh(self.casename)
        BCTab1.drawBCVizBoxes()
        
        # Create actors based on grid parts
        self.vtkActorList = []
        vtkgrids = self.data.buildVTKGrids()
        itag = 0
        for vtkgrid in vtkgrids:
            mapper = vtk.vtkDataSetMapper()
            mapper.SetInput(vtkgrid)
            self.vtkActorList.append(vtk.vtkActor())
            self.vtkActorList[-1].SetMapper(mapper)
            # set visibility of volumes to false on startup
            if int(self.data.elemGroups[itag]) < 0:
                self.vtkActorList[-1].VisibilityOff()
            itag = itag + 1
            
        #visualize the grid
        for actor in self.vtkActorList:
            self.ren.AddActor(actor)

        self.show()
        self.iren.Initialize()

    def addItems(self, parent):
        column = 0
        clients_item = self.addParent(parent, column, 'Clients', 'data Clients')
        vendors_item = self.addParent(parent, column, 'Vendors', 'data Vendors')
        time_period_item = self.addParent(parent, column, 'Time Period', 'data Time Period')

        
        self.addChild(clients_item, column, 'Type A', 'data Type A')
        self.addChild(clients_item, column, 'Type B', 'data Type B')

        self.addChild(vendors_item, column, 'Mary', 'data Mary')
        self.addChild(vendors_item, column, 'Arnold', 'data Arnold')

        self.addChild(time_period_item, column, 'Init', 'data Init')
        self.addChild(time_period_item, column, 'End', 'data End')

    def addParent(self, parent, column, title, data):
        item = QtGui.QTreeWidgetItem(parent, [title])
        item.setData(column, QtCore.Qt.UserRole, data)
        item.setChildIndicatorPolicy(QtGui.QTreeWidgetItem.ShowIndicator)
        item.setExpanded (True)
        return item

    def addChild(self, parent, column, title, data):
        item = QtGui.QTreeWidgetItem(parent, [title])
        item.setData(column, QtCore.Qt.UserRole, data)
        item.setCheckState (column, QtCore.Qt.Unchecked)
        return item

    def handleChanged(self, item, column):
        if item.checkState(column) == QtCore.Qt.Checked:
            print "checked", item, item.text(column)
        if item.checkState(column) == QtCore.Qt.Unchecked:
            print "unchecked", item, item.text(column)
    
    def visibilityChanged(self, vtkGridPartId, factagId, state):
        if state == QtCore.Qt.Checked:
            self.vtkActorList[vtkGridPartId].VisibilityOn()
            self.iren.Render()
            print "visibility On for factag: ", factagId
        else:
            self.vtkActorList[vtkGridPartId].VisibilityOff()
            self.iren.Render()
            print "visibility Off for factag: ", factagId

    def save(self):
        filename = "proteusCFD.state"
        pickle.dump(self.state, open(filename, "wb"))
        
    def load(self):
        filename = "proteusCFD.state"
        self.state = pickle.load(open(filename, "rb"))
        print self.state

    def selectFile(self):
        meshFile = QtGui.QFileDialog.getOpenFileName()
        parts = meshFile.split('.')
        self.casename = parts[0]
        print 'Selected casename ' + self.casename

        
        
if __name__ == "__main__":
 
    app = QtGui.QApplication(sys.argv)
 
    window = MainWindow()
 
    sys.exit(app.exec_())
