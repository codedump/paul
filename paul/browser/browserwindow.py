from PyQt4 import QtGui, QtCore

from paul.browser.matplotlibwidget import MatplotlibWidget
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar

import logging
log = logging.getLogger (__name__)

class BrowserWindow (QtGui.QMainWindow):
    def __init__ (self, start_path="/home/florin"):
        QtGui.QMainWindow.__init__(self)
        self.setWindowTitle ("Paul Browser")

        # some defaults
        self.start_path = start_path

        # create UI elements
        self.uiInit()

        # model for the file system (dir tree)
        self.fileSys = QtGui.QFileSystemModel()
        self.fileSys.setRootPath(QtCore.QDir.currentPath())
        self.fileSys.setFilter (QtCore.QDir.AllDirs | QtCore.QDir.Dirs | QtCore.QDir.Files | QtCore.QDir.NoDotAndDotDot)
        self.fileSys.setNameFilters ("*.ibw")
        self.fileSys.setNameFilterDisables (False)
        self.fileTree.setModel (self.fileSys)
        self.fileTree.activated.connect(self.dirActivated)
        

        # model for ibw files in the 2nd list box
        self.waveList.setRootIndex (self.fileSys.index(QtCore.QDir.currentPath()))
        self.waveList.setModel (self.fileSys)
        self.waveList.activated.connect(self.waveActivated)
        self.waveList.selectionModel().selectionChanged.connect(self.waveSelected)

        self.dirActivated (self.fileSys.index(self.start_path))

    def uiInit(self):
        self.main_frame = QtGui.QWidget()
        self.setCentralWidget (self.main_frame)

        self.hbox = QtGui.QHBoxLayout()
        self.main_frame.setLayout (self.hbox)

        #self.splitter = QtGui.QSplitter(self.main_frame)
        #self.splitter.setOrientation(QtCore.Qt.Horizontal)
        #self.splitter.setObjectName("splitter")

        self.fileTree = QtGui.QTreeView()
        self.waveList = QtGui.QListView()

        self.vbox_plot = QtGui.QVBoxLayout()
        self.plotCanvas = MatplotlibWidget()
        self.plotTools = NavigationToolbar(self.plotCanvas, self)
        
        for w in [ self.plotTools, self.plotCanvas ]:
            self.vbox_plot.addWidget (w)

        #self.plotDock = QtGui.QDockWidget(self.splitter)
        #self.plotDock.setEnabled(True)
        #self.plotDock.setFeatures(QtGui.QDockWidget.DockWidgetFloatable)

        for w in [ self.fileTree, self.waveList ]:
            self.hbox.addWidget (w)
            #self.hbox.setAlignment (w, QtCore.Qt.AlignVCenter)
        self.hbox.addLayout (self.vbox_plot)

    def ui_init_lists(self):
        pass

    def ui_init_plots(self):
        pass    


    # called when user selected an entry from the dir list
    @QtCore.pyqtSlot('QModelIndex')
    def dirActivated (self, index):
        self.waveList.setRootIndex (index)
        self.fileTree.scrollTo (index)
        self.fileTree.resizeColumnToContents(0)

    # called when an item in the waveList is activated (double-clicked)
    @QtCore.pyqtSlot('QModelIndex')
    def waveActivated (self, index):
        finfo = self.fileSys.fileInfo(index)
        fpath = self.fileSys.filePath(index)
        if finfo.isDir():
            self.waveList.setRootIndex (index)
            self.fileTree.scrollTo (index)
            self.fileTree.resizeColumnToContents(0)
        if finfo.isFile() and finfo.isReadable():
            log.info ("Loading %s" % fpath)
            self.plotCanvas.plotFile (fpath)

    # called when the selection changed in the waveList
    @QtCore.pyqtSlot ('QItemSelection', 'QItemSelection')
    def waveSelected (self, sel, unsel):
        if (sel.isEmpty()):
            return
        ind = sel.indexes()[0]
        log.debug ("Selected %s" % self.fileSys.filePath(ind))
        if self.fileSys.fileInfo(ind).isFile():
            log.debug ("Passing on %s" % self.fileSys.filePath(ind))
            self.waveActivated (ind)
