from PyQt4 import QtGui, QtCore

from paul.viewer.viewerwindow import ViewerWindow
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
        self.initViewer()

        # move to default start directory
        self.dirActivated (self.filesys.index(self.start_path))


    def initMainFrame(self):
        self.hbox = QtGui.QHBoxLayout()
        self.main_frame = QtGui.QWidget()
        self.main_frame.setLayout (self.hbox)
        self.setCentralWidget (self.main_frame)
        self.splitter = QtGui.QSplitter(self.main_frame)
        self.hbox.addWidget (self.splitter)


    def initBrowser(self):
        # model for the file system (dir tree)
        self.filesys = QtGui.QFileSystemModel()
        self.filesys.setRootPath(QtCore.QDir.currentPath())
        self.filesys.setFilter (QtCore.QDir.AllDirs | QtCore.QDir.Dirs | QtCore.QDir.Files | QtCore.QDir.NoDotAndDotDot)
        self.filesys.setNameFilters ("*.ibw")
        self.filesys.setNameFilterDisables (False)

        # the file system browser tree
        self.file_tree = QtGui.QTreeView(self.splitter)
        self.file_tree.activated.connect(self.dirActivated)
        self.file_tree.setModel (self.filesys)

        # list box for data files (waves)
        self.wave_list = QtGui.QListView(self.splitter)
        self.wave_list.setRootIndex (self.filesys.index(QtCore.QDir.currentPath()))
        self.wave_list.setModel (self.filesys)
        self.wave_list.activated.connect(self.waveActivated)
        self.wave_list.selectionModel().selectionChanged.connect(self.waveSelected)


    def initViewer(self):
        self.viewer = ViewerWindow(self.splitter)
        self.viewer.show()
        


    @QtCore.pyqtSlot('QModelIndex')
    def dirActivated (self, index):
        '''
        Called when user selected an entry from the dir list. The corresponding
        directory will be shown in the wave list.
        '''
        self.wave_list.setRootIndex (index)
        self.file_tree.scrollTo (index)
        self.file_tree.resizeColumnToContents(0)

    @QtCore.pyqtSlot('QModelIndex')
    def waveActivated (self, index):
        '''
        Called when an item in the waveList is activated (double-clicked).
        If the item is a wave, it is loaded and plotted. If it's a directory,
        it activated (i.e. descendend into).
        '''
        finfo = self.filesys.fileInfo(index)
        fpath = self.filesys.filePath(index)
        if finfo.isDir():
            self.wave_list.setRootIndex (index)
            self.file_tree.scrollTo (index)
            self.file_tree.resizeColumnToContents(0)
            self.last_path = str(fpath)
        if finfo.isFile() and finfo.isReadable():
            log.info ("Loading %s" % fpath)
            self.viewer.plotFile (fpath)

    # called when the selection changed in the waveList
    @QtCore.pyqtSlot ('QItemSelection', 'QItemSelection')
    def waveSelected (self, sel, unsel):
        '''
        Called when user selects (marks or unmarks) an item in the wave list.
        The general idea is to load and plot the first selected wave (if
        it's 2D), or all selected waves (if they are 1D).
        '''
        if (sel.isEmpty()):
            return
        ind = sel.indexes()[0]
        log.debug ("Selected %s" % self.filesys.filePath(ind))
        if self.filesys.fileInfo(ind).isFile():
            log.debug ("Passing on %s" % self.filesys.filePath(ind))
            self.waveActivated (ind)
