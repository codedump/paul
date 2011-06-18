# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '../../ui/browserwindow.ui'
#
# Created: Sa Jun 18 16:07:22 2011
#      by: The PyQt User Interface Compiler (pyuic) 3.18.1
#
# WARNING! All changes made in this file will be lost!


from qt import *
from paul.browser.mplwidget import *


class BrowserWindow(QMainWindow):
    def __init__(self,parent = None,name = None,fl = 0):
        QMainWindow.__init__(self,parent,name,fl)
        self.statusBar()

        if not name:
            self.setName("BrowserWindow")

        self.setSizePolicy(QSizePolicy(QSizePolicy.Preferred,QSizePolicy.Preferred,0,0,self.sizePolicy().hasHeightForWidth()))

        self.setCentralWidget(QWidget(self,"qt_central_widget"))

        self.Plotspace = MatplotlibWidget(self.centralWidget(),"Plotspace")
        self.Plotspace.setGeometry(QRect(10,10,560,420))
        self.Plotspace.setSizePolicy(QSizePolicy(QSizePolicy.Maximum,QSizePolicy.Maximum,0,0,self.Plotspace.sizePolicy().hasHeightForWidth()))
        self.Plotspace.setMinimumSize(QSize(50,50))



        self.languageChange()

        self.resize(QSize(580,453).expandedTo(self.minimumSizeHint()))
        self.clearWState(Qt.WState_Polished)


    def languageChange(self):
        self.setCaption(self.__tr("Paul Browser"))


    def __tr(self,s,c = None):
        return qApp.translate("BrowserWindow",s,c)
