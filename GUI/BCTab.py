from guiData import *
from PyQt4 import QtCore, QtGui
import functools

class BCTab(QtGui.QWidget):
    """This class represents the bc tab"""

    # pass in the gui data object and a function to change visibility of vtk parts
    def __init__(self, data, visibilityChangeFunc):
        self.data = data
        self.visibilityChangeFunc = visibilityChangeFunc
        super(BCTab, self).__init__()
        addBCButton = QtGui.QPushButton('Add Boundary Condition', self)
        addBCButton.clicked.connect(functools.partial(self.BoundaryConditionPopup, self.data))
        addBCButton.setToolTip('Add a new boundary condition')
        addBCButton.resize(addBCButton.sizeHint())
        addBCButton.move(10, 10)

    def clearBCVizBoxes(self):
        checkBoxes = self.findChildren(QtGui.QCheckBox)
        for i in reversed(range(0, len(checkBoxes))): 
            checkBoxes[i].setParent(None)
        
    def drawBCVizBoxes(self):
        # Clear the old boxes if they exist
        self.clearBCVizBoxes()
        # This loop assumes that the vtk parts are created in the same order as the factags
        irow = 0
        for tag in self.data.elemGroups:
            if int(tag) >= 0:
                checkViz = QtGui.QCheckBox('factag: ' + str(tag), self)
                checkViz.setChecked(True)
            else:
                checkViz = QtGui.QCheckBox('voltag: ' + str(tag), self)
                checkViz.setChecked(False) # hide volumes by default
            checkViz.stateChanged.connect(functools.partial(self.visibilityChangeFunc, int(irow), int(tag)))
            checkViz.move(20, 40+25*irow)
            irow = irow + 1
        self.bcwindow = None

    def BoundaryConditionPopup(self, data):
        print "Opening a new BC popup window..."
        self.bcwindow = BCPopup(data)
        self.bcwindow.setGeometry(QtCore.QRect(100, 100, 600, 300))
        self.bcwindow.show()

class BCPopup(QtGui.QWidget):
    """This class represents the bc popup for BC editing"""
    
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
        self.addBCButton.clicked.connect(self.creatBCObj)
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

    def createBCObj(self):
        print 'Apply clicked'
