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
        self.initMainFrame()
        self.initBrowser()

        # move to default start directory
        self.dirActivated (self.filesys.index(self.start_path))


    def initMainFrame(self):
        self.main_frame = QtGui.QWidget()
        self.setCentralWidget (self.main_frame)

        self.hbox = QtGui.QHBoxLayout()
        self.main_frame.setLayout (self.hbox)
        self.splitter = QtGui.QSplitter(self.main_frame)
        self.hbox.addWidget (self.splitter)

        self.file_tree = QtGui.QTreeView(self.splitter)
        self.wave_list = QtGui.QListView(self.splitter)

        self.vbox_plot = QtGui.QVBoxLayout()
        self.plot_canvas = MatplotlibWidget()
        self.plot_tools = NavigationToolbar(self.plot_canvas, self)
        
        for w in [ self.plot_tools, self.plot_canvas ]:
            self.vbox_plot.addWidget (w)

        #self.plotDock = QtGui.QDockWidget(self.splitter)
        #self.plotDock.setEnabled(True)
        #self.plotDock.setFeatures(QtGui.QDockWidget.DockWidgetFloatable)

        self.hbox.addLayout (self.vbox_plot)


    def initBrowser(self):
        # model for the file system (dir tree)
        self.filesys = QtGui.QFileSystemModel()
        self.filesys.setRootPath(QtCore.QDir.currentPath())
        self.filesys.setFilter (QtCore.QDir.AllDirs | QtCore.QDir.Dirs | QtCore.QDir.Files | QtCore.QDir.NoDotAndDotDot)
        self.filesys.setNameFilters ("*.ibw")
        self.filesys.setNameFilterDisables (False)
        self.file_tree.setModel (self.filesys)
        self.file_tree.activated.connect(self.dirActivated)        

        # model for ibw files in the 2nd list box
        self.wave_list.setRootIndex (self.filesys.index(QtCore.QDir.currentPath()))
        self.wave_list.setModel (self.filesys)
        self.wave_list.activated.connect(self.waveActivated)
        self.wave_list.selectionModel().selectionChanged.connect(self.waveSelected)

    def ui_init_plots(self):
        pass    


    # called when user selected an entry from the dir list
    @QtCore.pyqtSlot('QModelIndex')
    def dirActivated (self, index):
        self.wave_list.setRootIndex (index)
        self.file_tree.scrollTo (index)
        self.file_tree.resizeColumnToContents(0)

    # called when an item in the waveList is activated (double-clicked)
    @QtCore.pyqtSlot('QModelIndex')
    def waveActivated (self, index):
        finfo = self.filesys.fileInfo(index)
        fpath = self.filesys.filePath(index)
        if finfo.isDir():
            self.wave_list.setRootIndex (index)
            self.file_tree.scrollTo (index)
            self.file_tree.resizeColumnToContents(0)
        if finfo.isFile() and finfo.isReadable():
            log.info ("Loading %s" % fpath)
            self.plot_canvas.plotFile (fpath)

    # called when the selection changed in the waveList
    @QtCore.pyqtSlot ('QItemSelection', 'QItemSelection')
    def waveSelected (self, sel, unsel):
        if (sel.isEmpty()):
            return
        ind = sel.indexes()[0]
        log.debug ("Selected %s" % self.filesys.filePath(ind))
        if self.filesys.fileInfo(ind).isFile():
            log.debug ("Passing on %s" % self.filesys.filePath(ind))
            self.waveActivated (ind)
