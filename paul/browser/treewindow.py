from PyQt4 import QtGui, QtCore

from paul.viewer.viewerwindow import ViewerWindow
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar

import os.path
import logging
log = logging.getLogger (__name__)

class TreeWindow (QtGui.QMainWindow):

    # Signal emited when user selected one or more files for plotting.
    wavesSelected = QtCore.pyqtSignal('QStringList')

    current_path = ''

    def __init__ (self, start_path="/home/florin"):
        QtGui.QMainWindow.__init__(self)
        self.setWindowTitle ("Paul Browser")

        # create UI elements
        self.initMainFrame()
        self.initBrowser()

        # move to default start directory
        self.setRoot (start_path)


    def initMainFrame(self):
        '''Initializes the inner frame of the tree-view window (main browser feature)'''
        self.vbox = QtGui.QVBoxLayout()
        self.main_frame = QtGui.QWidget()
        self.main_frame.setLayout (self.vbox)
        self.setCentralWidget (self.main_frame)


    def initBrowser(self):
        '''Initializes browser specific stuff (tree-view and buttons)'''

        # model for the file system (dir tree)
        self.filesys = QtGui.QFileSystemModel()
        self.filesys.setRootPath(QtCore.QDir.currentPath())
        self.filesys.setFilter (QtCore.QDir.AllDirs | QtCore.QDir.Dirs | QtCore.QDir.Files | QtCore.QDir.NoDotAndDotDot)
        self.filesys.setNameFilters ("*.ibw")
        self.filesys.setNameFilterDisables (False)

        # create general tree browser layout:
        # one giant VBox, holding the buttons (toolbar?) in the upper part
        # and the file tree in the lower part.
        self.tools_hbox = QtGui.QHBoxLayout()
        self.vbox.addLayout (self.tools_hbox)
        self.tools_btnUp = QtGui.QPushButton ("&Up")
        self.tools_btnUp.show()
        self.tools_btnUp.setSizePolicy (QtGui.QSizePolicy.Minimum,
                                        QtGui.QSizePolicy.Minimum)
        self.tools_btnUp.clicked.connect (self.treeUp)
        self.tools_hbox.addWidget (self.tools_btnUp)

        # the file system browser tree
        self.filetree = QtGui.QTreeView()
        self.vbox.addWidget (self.filetree)
        self.filetree.activated.connect(self.itemActivated)
        self.filetree.setModel (self.filesys)
        self.filetree.selectionModel().selectionChanged.connect(self.itemSelected)

        
    def setRoot (self, path):
        '''
        Sets the specified path as the current root path. Typically
        called when the user double-clicked an entry in the tree view.
        '''
        self.current_path = os.path.abspath(path)
        self.statusBar().showMessage (self.current_path)
        index = self.filesys.index (self.current_path)
        self.filetree.setRootIndex (index)
        self.filetree.scrollTo (index)

        
    @QtCore.pyqtSlot()
    def treeUp(self):
        '''
        Moves the root of the tree view one item up (i.e.
        sets the parent of the current root as the new root).
        '''
        new_path = os.path.dirname(self.current_path)
        log.debug ("Moving to parent: %s" % new_path)
        self.setRoot (new_path)


    @QtCore.pyqtSlot('QModelIndex')
    def waveActivated (self, index):
        '''
        Called when an item is activated and that item is a
        regular file (i.e. probably a wave).
        If the item is a wave, it is loaded and plotted.
        '''
        finfo = self.filesys.fileInfo(index)
        fpath = self.filesys.filePath(index)
        log.debug ("Activated file: %s" % fpath)


    @QtCore.pyqtSlot('QModelIndex')
    def itemActivated (self, index):
        '''
        Called when an item in the item tree is activated (double-clicked).
        If the item is a wave, it is loaded and plotted. If it's a directory,
        it activated (i.e. descendend into).
        '''
        finfo = self.filesys.fileInfo(index)
        if finfo.isDir():
            self.setRoot (str(self.filesys.filePath(index)))
        if finfo.isFile() and finfo.isReadable():
            self.waveActivated (index)


    # called when the selection changed in the waveList
    @QtCore.pyqtSlot ('QItemSelection', 'QItemSelection')
    def itemSelected (self, sel, unsel):
        '''
        Called when user selects (marks or unmarks) an item in the wave list.
        The general idea is to load and plot the first selected wave (if
        it's 2D), or all selected waves (if they are 1D).
        '''
        if (sel.isEmpty()):
            return

        ind = sel.indexes()[0]
        finfo = self.filesys.fileInfo(ind)
        fpath = self.filesys.filePath(ind)

        if finfo.isFile():
            log.debug ("Selected file: %s" % str(fpath))
            self.wavesSelected.emit([str(fpath)])
        elif finfo.isDir():
            log.debug ("Selected dir: %s" % str(fpath))
        else:
            log.debug ("...what to do with %s?" % str(fpath))
