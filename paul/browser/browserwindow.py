from PyQt4 import QtGui, QtCore
from paul.browser.ui_browserwindow import Ui_BrowserWindow

class BrowserWindow (QtGui.QMainWindow, Ui_BrowserWindow):
    def __init__ (self):
        QtGui.QMainWindow.__init__(self)
        Ui_BrowserWindow.__init__(self)
        self.setupUi (self)
        
        self.fileSys = QtGui.QFileSystemModel()
        self.fileSys.setRootPath(QtCore.QDir.currentPath())
        
        self.fileTree.setModel (self.fileSys)
