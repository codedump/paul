from PyQt4 import QtGui, QtCore
from paul.browser.ui_browserwindow import Ui_BrowserWindow

class BrowserWindow (QtGui.QMainWindow, Ui_BrowserWindow):
    def __init__ (self):
        QtGui.QMainWindow.__init__(self)
        Ui_BrowserWindow.__init__(self)

        self.setupUi (self)

        # model for the file system (dir tree)
        self.fileSys = QtGui.QFileSystemModel()
        self.fileSys.setRootPath("~/")

        # model for ibw files in the 2nd list box
        self.fileIbw = QtGui.QFileSystemModel()
        self.fileIbw.setFilter (QtCore.QDir.Files)
        self.fileIbw.setNameFilters ("*.ibw")

        self.fileTree.setModel (self.fileSys)
        self.fileList.setModel (self.fileIbw)
        self.fileTree.setColumnWidth (0, 200)

        self.fileTree.activated.connect(self.dirSelected)

    # called when user selected an entry from the dir list
    @QtCore.pyqtSlot('QModelIndex')
    def dirSelected (self, index):
        print "selected: %s" % self.fileSys.filePath(index)
        self.fileIbw.setRootPath (self.fileSys.filePath(index))
        
