#!/usr/bin/python

from PyQt4 import QtCore

import logging
log = logging.getLogger (__name__)

class WaveModel (QtCore.QAbstractItemModel):

    waveList = ['one', 'two', 'three']

    def __init__ (self, parent = 0):
        QtCore.QObject.__init__(self)
        self.basePath = QtCore.QString("~")

    def columnCount (self, parent):
        return 1

    def rowCount (self, parent):
        return len(self.waveList)

    def index(self, row, col, ind):
        ind = self.createIndex(row, col)
        return ind

    def data (self, index, role):
        if (role == QtCore.Qt.DisplayRole):
            return self.waveList[index.row()]

    @QtCore.pyqtSlot('QString')
    def setBasePath (self, np):
        log.debug ("Selected: %s" % np)
        self.basePath = np
