# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '../../ui/browserwindow.ui'
#
# Created: Thu May 10 14:59:09 2012
#      by: PyQt4 UI code generator 4.7.3
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

class Ui_BrowserWindow(object):
    def setupUi(self, BrowserWindow):
        BrowserWindow.setObjectName("BrowserWindow")
        BrowserWindow.resize(444, 290)
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
        self.horizontalLayout = QtGui.QHBoxLayout(self.widget)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.splitter = QtGui.QSplitter(self.widget)
        self.splitter.setOrientation(QtCore.Qt.Horizontal)
        self.splitter.setObjectName("splitter")
        self.fileTree = QtGui.QTreeView(self.splitter)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.fileTree.sizePolicy().hasHeightForWidth())
        self.fileTree.setSizePolicy(sizePolicy)
        self.fileTree.setObjectName("fileTree")
        self.waveList = QtGui.QListView(self.splitter)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Preferred, QtGui.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.waveList.sizePolicy().hasHeightForWidth())
        self.waveList.setSizePolicy(sizePolicy)
        self.waveList.setObjectName("waveList")
        self.plotDock = QtGui.QDockWidget(self.splitter)
        self.plotDock.setEnabled(True)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.plotDock.sizePolicy().hasHeightForWidth())
        self.plotDock.setSizePolicy(sizePolicy)
        self.plotDock.setMinimumSize(QtCore.QSize(118, 272))
        self.plotDock.setBaseSize(QtCore.QSize(300, 300))
        self.plotDock.setFeatures(QtGui.QDockWidget.DockWidgetFloatable)
        self.plotDock.setObjectName("plotDock")
        self.plotCanvas = MatplotlibWidget()
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.plotCanvas.sizePolicy().hasHeightForWidth())
        self.plotCanvas.setSizePolicy(sizePolicy)
        self.plotCanvas.setObjectName("plotCanvas")
        self.plotDock.setWidget(self.plotCanvas)
        self.horizontalLayout.addWidget(self.splitter)
        BrowserWindow.setCentralWidget(self.widget)

        self.retranslateUi(BrowserWindow)
        QtCore.QMetaObject.connectSlotsByName(BrowserWindow)

    def retranslateUi(self, BrowserWindow):
        BrowserWindow.setWindowTitle(QtGui.QApplication.translate("BrowserWindow", "Paul Browser", None, QtGui.QApplication.UnicodeUTF8))

from matplotlibwidget import MatplotlibWidget
