#!/usr/bin/env python
 
import sys
import pickle # we use pickle to save run information, etc.
import vtk
from guiData import *
from PyQt4 import QtCore, QtGui
from vtk.qt4.QVTKRenderWindowInteractor import QVTKRenderWindowInteractor
from defines import *
import functools


class MainWindow(QtGui.QMainWindow):
 
    def __init__(self, parent = None):
        QtGui.QMainWindow.__init__(self, parent)

        self.state = {"test", "stuff", "other"}
        self.data = GuiData()
        #read the mesh
        self.data.loadMesh("bump")
         
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
        BCTab = QtGui.QWidget()
        self.tabs.addTab(BCTab, "Boundary Conditions")
        addBCButton = QtGui.QPushButton('Add Boundary Condition', BCTab)
        addBCButton.clicked.connect(self.BoundaryConditionPopup)
        addBCButton.setToolTip('Add a new boundary condition')
        addBCButton.resize(addBCButton.sizeHint())
        addBCButton.move(10, 10)
        irow = 0
        for tag in self.data.elemGroups:
            if int(tag) >= 0:
                checkViz = QtGui.QCheckBox('factag: ' + str(tag), BCTab)
                checkViz.setChecked(True)
            else:
                checkViz = QtGui.QCheckBox('voltag: ' + str(tag), BCTab)
                checkViz.setChecked(False) # hide volumes by default
            checkViz.stateChanged.connect(functools.partial(self.visibilityChanged, int(irow), int(tag)))
            checkViz.move(20, 40+25*irow)
            irow = irow + 1
        self.bcwindow = None
        
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
        
        self.menubar = self.menuBar()
        fileMenu = self.menubar.addMenu('&File')
        fileMenu.addAction(exitAction)
        fileMenu.addAction(saveAction)
        fileMenu.addAction(loadAction)
        
        self.statusBar().showMessage('Waiting...')
        
        #Setup our vtk widget
        self.vtkWidget = QVTKRenderWindowInteractor(self.frame)
        self.vl.addWidget(self.vtkWidget, 0, 300)

        self.ren = vtk.vtkRenderer()
        self.ren.SetBackground(.5,.8,.5)       
        self.ren.ResetCamera()

        self.vtkWidget.GetRenderWindow().AddRenderer(self.ren)
        self.iren = self.vtkWidget.GetRenderWindow().GetInteractor()

        # Create actors based on grid parts
        self.vtkActorList = []
        vtkgrids = self.data.buildVTKGrids()
        itag = 0
        for vtkgrid in vtkgrids:
            mapper = vtk.vtkDataSetMapper()
            mapper.SetInput(vtkgrid)
            self.vtkActorList.append(vtk.vtkActor())
            self.vtkActorList[-1].SetMapper(mapper)
            if int(self.data.elemGroups[itag]) < 0:
                self.vtkActorList[-1].VisibilityOff()
            
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

    def BoundaryConditionPopup(self):
        print "Opening a new BC popup window..."
        self.bcwindow = BCPopup(self.data)
        self.bcwindow.setGeometry(QtCore.QRect(100, 100, 600, 300))
        self.bcwindow.show()

class BCPopup(QtGui.QWidget):
    def __init__(self, data):
        QtGui.QWidget.__init__(self)
        self.numBCsSet = len(data.elemGroups)
        self.bcIdSelected = 0
        self.bcTypeSelected = 'NULL'

        # setup layout
        self.vl = QtGui.QGridLayout()
        self.vl.setSpacing(10)
        self.setLayout(self.vl)

        # List of factags (i.e. BCs to set)
        self.bcList = QtGui.QComboBox()
        for tag in data.elemGroups:
            self.bcList.addItem('factag: ' + str(tag))
        self.bcList.activated[str].connect(self.bcListActivate)
        self.vl.addWidget(self.bcList,0,0)
        
        # List of bc types available
        self.typeList = QtGui.QComboBox()
        for type in data.bcTypes:
            self.typeList.addItem(type)
        self.typeList.activated[str].connect(self.bcTypeActivate)
        self.vl.addWidget(self.typeList,1,0)

        # List of current BCs set
        bcCols = 4
        bcRows = self.numBCsSet
        self.bcTable = QtGui.QTableWidget()
        self.bcTable.setRowCount(self.numBCsSet)
        self.bcTable.setColumnCount(bcCols)
        self.bcTable.resize(self.bcTable.sizeHint())
        irow = 0
        for tag in data.elemGroups:
            self.bcTable.setItem(irow,1, QtGui.QTableWidgetItem(str(tag)))
            irow = irow + 1
            
        self.bcTable.setHorizontalHeaderLabels(("BoundaryID; Factag; head3; Nickname").split(";"))
        self.bcTable.show()
        self.vl.addWidget(self.bcTable,2,1)


        self.addBCButton = QtGui.QPushButton('Apply')
        #addBCButton.clicked.connect()
        self.addBCButton.setToolTip('Apply Boundary Condition')
        self.addBCButton.resize(self.addBCButton.sizeHint())
        self.addBCButton.setStyleSheet("background-color: #5BC85B")
        self.vl.addWidget(self.addBCButton,3,1)

    def bcListActivate(self, text):
        print 'Combo selection changed to ' + text
        self.bcIdSelected = text

    def bcTypeActivate(self, text):
        print 'Combo selection changed to ' + text
        self.bcTypeSelected = text

    def treeChanged(self, item, column):
        print 'tree changed'
        
        
if __name__ == "__main__":
 
    app = QtGui.QApplication(sys.argv)
 
    window = MainWindow()
 
    sys.exit(app.exec_())
