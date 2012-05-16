from PyQt4 import QtGui, QtCore

from paul.browser.ui_browserwindow import Ui_BrowserWindow

import logging
log = logging.getLogger (__name__)

class BrowserWindow (QtGui.QMainWindow, Ui_BrowserWindow):
    def __init__ (self, start_path="/home/florin"):
        QtGui.QMainWindow.__init__(self)
        Ui_BrowserWindow.__init__(self)

        # some defaults
        self.start_path = start_path

        self.setupUi (self)

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
