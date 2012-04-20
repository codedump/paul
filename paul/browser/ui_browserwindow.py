# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '../../ui/browserwindow.ui'
#
# Created: Fri Apr 20 11:03:40 2012
#      by: PyQt4 UI code generator 4.7.3
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

class Ui_BrowserWindow(object):
    def setupUi(self, BrowserWindow):
        BrowserWindow.setObjectName("BrowserWindow")
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
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
        self.gridLayout = QtGui.QGridLayout(self.widget)
        self.gridLayout.setObjectName("gridLayout")
        self.splitter = QtGui.QSplitter(self.widget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.splitter.sizePolicy().hasHeightForWidth())
        self.splitter.setSizePolicy(sizePolicy)
        self.splitter.setOrientation(QtCore.Qt.Horizontal)
        self.splitter.setObjectName("splitter")
        self.fileTree = QtGui.QTreeView(self.splitter)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.fileTree.sizePolicy().hasHeightForWidth())
        self.fileTree.setSizePolicy(sizePolicy)
        self.fileTree.setBaseSize(QtCore.QSize(100, 100))
        self.fileTree.setObjectName("fileTree")
        self.fileList = QtGui.QListView(self.splitter)
        self.fileList.setEnabled(True)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.fileList.sizePolicy().hasHeightForWidth())
        self.fileList.setSizePolicy(sizePolicy)
        self.fileList.setBaseSize(QtCore.QSize(100, 100))
        self.fileList.setObjectName("fileList")
        self.plotDock = QtGui.QDockWidget(self.splitter)
        self.plotDock.setEnabled(True)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.plotDock.sizePolicy().hasHeightForWidth())
        self.plotDock.setSizePolicy(sizePolicy)
        self.plotDock.setMinimumSize(QtCore.QSize(70, 38))
        self.plotDock.setObjectName("plotDock")
        self.plotDockContents = QtGui.QWidget()
        self.plotDockContents.setObjectName("plotDockContents")
        self.gridLayout_2 = QtGui.QGridLayout(self.plotDockContents)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.plotCanvas = MatplotlibWidget(self.plotDockContents)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.plotCanvas.sizePolicy().hasHeightForWidth())
        self.plotCanvas.setSizePolicy(sizePolicy)
        self.plotCanvas.setMinimumSize(QtCore.QSize(0, 0))
        self.plotCanvas.setObjectName("plotCanvas")
        self.gridLayout_2.addWidget(self.plotCanvas, 0, 0, 1, 1)
        self.plotDock.setWidget(self.plotDockContents)
        self.gridLayout.addWidget(self.splitter, 0, 0, 1, 1)
        BrowserWindow.setCentralWidget(self.widget)

        self.retranslateUi(BrowserWindow)
        QtCore.QMetaObject.connectSlotsByName(BrowserWindow)

    def retranslateUi(self, BrowserWindow):
        BrowserWindow.setWindowTitle(QtGui.QApplication.translate("BrowserWindow", "Paul Browser", None, QtGui.QApplication.UnicodeUTF8))

from matplotlibwidget import MatplotlibWidget
