from PyQt4 import QtGui, QtCore

from paul.browser.ui_browserwindow import Ui_BrowserWindow
from paul.browser.wavemodel import WaveModel

import logging
log = logging.getLogger (__name__)

class BrowserWindow (QtGui.QMainWindow, Ui_BrowserWindow):
    def __init__ (self):
        QtGui.QMainWindow.__init__(self)
        Ui_BrowserWindow.__init__(self)

        self.setupUi (self)

        # model for the file system (dir tree)
        self.fileSys = QtGui.QFileSystemModel()
        self.fileSys.setRootPath("/home/florin")
        self.fileTree.setModel (self.fileSys)
        self.fileTree.activated.connect(self.dirSelected)

        # model for ibw files in the 2nd list box
        self.waveMod = WaveModel()
        self.waveList.setModel (self.waveMod)

    # called when user selected an entry from the dir list
    @QtCore.pyqtSlot('QModelIndex')
    def dirSelected (self, index):
        self.waveMod.setBasePath (self.fileSys.filePath(index))
        
