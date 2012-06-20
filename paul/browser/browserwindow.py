from PyQt4 import QtGui, QtCore

from paul.viewer.viewerwindow import ViewerWindow
from paul.browser.treewindow import TreeWindow
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar

import logging
log = logging.getLogger (__name__)

class BrowserWindow (QtGui.QMainWindow):
    def __init__ (self, start_path='~'):
        QtGui.QMainWindow.__init__(self)
        self.setWindowTitle ("Paul Browser")

        # create UI elements
        self.initMainFrame()
        self.initBrowser()
        self.initViewer()

        for w in [ self.viewer, self.tree ]:
            w.vbox.setContentsMargins (0, 0, 0, 0)

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
