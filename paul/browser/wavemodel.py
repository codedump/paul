#!/usr/bin/python

from PyQt4 import QtCore

import logging
log = logging.getLogger (__name__)

#
# OBSOLETE (for now...)
#

#
# This class implements a model for a list of (viewable) waves,
# i.e. data files that can be read and plotted. The idea
# is that the user, having selected a "base path" by other means
# (i.e. a tree view), gets to see all viewable items (from a Paul
# point of view) from that path component.
#
# The path component (basePath) can be a directory or a file.
# In the most general case, a directory, the directory contents
# (mostly files) following a general naming pattern, known to be
# readable by Paul (e.g. .IBW files) will be shown.
#
# If the path component is a file (i.e. an IgorPro .IBW, or Elmitec
# .DAT file), or a multi-wave construct (i.e. a HDF5-container), then
# the contents of that container will be displayed... (as soon as it's
# implemented :-) )
#
class WaveModel (QtCore.QAbstractItemModel):

    def __init__ (self, parent = 0):
        QtCore.QAbstractItemModel.__init__(self)
        self.basePath = QtCore.QString("")
        self.waveList = QtCore.QStringList()
        self.fmtList = QtCore.QStringList()

    def columnCount (self, parent):
        return 1

    def rowCount (self, parent):
        return self.waveList.count()

    def index(self, row, col, parent):
        # FIXME: Directly indexing a list member... *ugh* this is very ugly :-(
        #        There _has_ to be a better way, i.e. direct access using
        #        some magic of QModelIndex, but I haven't figured it out yet.
        if not parent.isValid():
            return self.createIndex (row, col, (str(self.waveList[row])))
        else:
            return QtCore.QModelIndex()

    def parent (self, child_index):
        # Return the parent of child_index.
        # As of now, only top-level items are present,
        # so we return an invalid index.
        return QtCore.QModelIndex()

    def data (self, index, role):
        if (role == QtCore.Qt.DisplayRole):
            return index.internalPointer()
        else:
            return None

    # Called when a new base path is to be reloaded.
    @QtCore.pyqtSlot('QString')
    def setBasePath (self, np):
        log.debug ("Selected: %s" % np)
        self.beginResetModel()
        self.basePath = np
        self.baseDir = QtCore.QDir (self.basePath)
        self.dirWatcher = QtCore.QFileSystemWatcher()
        self.dirWatcher.addPath (self.basePath)
        self.dirWatcher.directoryChanged.connect(self.setBasePath)
        self.waveList = self.baseDir.entryList (self.fmtList)
        log.debug ("%d items in %s" % (self.waveList.count(), self.baseDir.absolutePath()))
        self.endResetModel()
