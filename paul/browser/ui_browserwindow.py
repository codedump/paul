# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '../../ui/browserwindow.ui'
#
# Created: Thu Apr 19 15:16:24 2012
#      by: PyQt4 UI code generator 4.7.3
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

class Ui_BrowserWindow(object):
    def setupUi(self, BrowserWindow):
        BrowserWindow.setObjectName("BrowserWindow")
        BrowserWindow.resize(786, 402)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(BrowserWindow.sizePolicy().hasHeightForWidth())
        BrowserWindow.setSizePolicy(sizePolicy)
        BrowserWindow.setProperty("dockWindowsMovable", True)
        self.widget = QtGui.QWidget(BrowserWindow)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.widget.sizePolicy().hasHeightForWidth())
        self.widget.setSizePolicy(sizePolicy)
        self.widget.setObjectName("widget")
        self.splitter = QtGui.QSplitter(self.widget)
        self.splitter.setGeometry(QtCore.QRect(0, 20, 590, 350))
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.splitter.sizePolicy().hasHeightForWidth())
        self.splitter.setSizePolicy(sizePolicy)
        self.splitter.setMaximumSize(QtCore.QSize(590, 350))
        self.splitter.setAutoFillBackground(False)
        self.splitter.setFrameShadow(QtGui.QFrame.Plain)
        self.splitter.setOrientation(QtCore.Qt.Horizontal)
        self.splitter.setOpaqueResize(True)
        self.splitter.setObjectName("splitter")
        self.fileTree = QtGui.QTreeView(self.splitter)
        self.fileTree.setObjectName("fileTree")
        self.fileList = QtGui.QListView(self.splitter)
        self.fileList.setEnabled(True)
        self.fileList.setObjectName("fileList")
        self.plotDock = QtGui.QDockWidget(self.splitter)
        self.plotDock.setObjectName("plotDock")
        self.plotDockContents = QtGui.QWidget()
        self.plotDockContents.setObjectName("plotDockContents")
        self.plotCanvas = MatplotlibWidget(self.plotDockContents)
        self.plotCanvas.setGeometry(QtCore.QRect(9, 9, 52, 50))
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.plotCanvas.sizePolicy().hasHeightForWidth())
        self.plotCanvas.setSizePolicy(sizePolicy)
        self.plotCanvas.setMinimumSize(QtCore.QSize(50, 50))
        self.plotCanvas.setObjectName("plotCanvas")
        self.plotDock.setWidget(self.plotDockContents)
        BrowserWindow.setCentralWidget(self.widget)

        self.retranslateUi(BrowserWindow)
        QtCore.QMetaObject.connectSlotsByName(BrowserWindow)

    def retranslateUi(self, BrowserWindow):
        BrowserWindow.setWindowTitle(QtGui.QApplication.translate("BrowserWindow", "Paul Browser", None, QtGui.QApplication.UnicodeUTF8))

from ui_matplotlibwidget import MatplotlibWidget
