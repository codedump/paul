from PyQt4 import QtGui, QtCore

from paul.viewer.viewerwindow import ViewerWindow
from paul.browser.treewindow import TreeWindow
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
import paul.loader.igor as igor

import os.path

import logging
log = logging.getLogger (__name__)

class BrowserWindow (QtGui.QMainWindow):
    def __init__ (self, start_path='~'):
        QtGui.QMainWindow.__init__(self)
        self.setWindowTitle ("Paul Browser")
        self.resize (800, 500)

        # create UI elements
        self.initMainFrame()
        self.initBrowser()
        self.initViewer()

        # inter-component connections
        self.tree.wavesSelected.connect(self.viewer.plotFiles)

        # move to default start directory
        self.tree.setRoot (start_path)


    def initMainFrame(self):
        self.hbox = QtGui.QHBoxLayout()
        self.hbox.setSpacing (0)
        self.hbox.setContentsMargins(0, 0, 0, 0)
        self.main_frame = QtGui.QWidget()
        self.main_frame.setLayout (self.hbox)
        self.setCentralWidget (self.main_frame)
        self.splitter = QtGui.QSplitter(self.main_frame)
        self.hbox.addWidget (self.splitter)


    def initBrowser(self):
        # model for the file system (dir tree)
        self.tree = TreeWindow()
        self.tree.setParent (self.splitter)


    def initViewer(self):
        self.viewer = ViewerWindow()
        self.viewer.setParent (self.splitter)


    @QtCore.pyqtSlot('QStringList')
    def plotFiles(wfiles):
        # ...
        pass

    #
    # Some useful browser properties
    #

    @property
    def root(self):
        '''
        The root component of the browser.
        '''
        return self.tree.current_path

    @root.setter
    def root(self, path):
        '''
        Sets the root of the browser.
        '''
        return self.tree.setRoot (path)

    @property
    def _sel(self):
        '''
        Same as *BrowserWindow.sel*, only the selection returns full
        a list of full paths instead of wave names. This also works correctly
        for waves with similar names in different directories.
        '''
        return self.tree.selected_paths

    @property
    def sel(self):
        '''
        The currently selected waves as a dictionary of full file paths,
        with key indices being the basename of the file path.
        (Note that files with the same basename in different directories
        will misbehave. Use the *BrowserWindow._sel* property instead.)
        '''
        return dict({os.path.basename(i): os.path.abspath(i) for i in self._sel})

    @property
    def wav(self):
        '''
        A dictionary of Wave() objects, using file base names (i.e. 
        the part of the path after the last '/') as keys. It will misbehave
        if different waves in different directories have the same file name.
        For this purpose, use *BrowserWindow._wav* instead.
        '''
        return dict({i: igor.load(j) for i,j in self.sel})


    @property
    def _wav(self):
        '''
        The currently selected waves as a list of Wave() objects.
        Works correctly for similary named waves in different directories.
        '''
        return [igor.load(i) for i in self._sel]
