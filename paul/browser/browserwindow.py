# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '../../ui/browserwindow.ui'
#
# Created: Fri Jun 17 14:34:09 2011
#      by: The PyQt User Interface Compiler (pyuic) 3.17.4
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


        self.setCentralWidget(QWidget(self,"qt_central_widget"))

        self.Dirlist = QListView(self.centralWidget(),"Dirlist")
        self.Dirlist.addColumn(self.__tr("Column 1"))
        self.Dirlist.setGeometry(QRect(11,11,150,420))

        self.Plotspace = MatplotlibWidget(self.centralWidget(),"Plotspace")
        self.Plotspace.setGeometry(QRect(169,11,400,418))
        self.Plotspace.setSizePolicy(QSizePolicy(QSizePolicy.Preferred,QSizePolicy.Preferred,0,0,self.Plotspace.sizePolicy().hasHeightForWidth()))
        self.Plotspace.setMinimumSize(QSize(50,50))



        self.languageChange()

        self.resize(QSize(580,453).expandedTo(self.minimumSizeHint()))
        self.clearWState(Qt.WState_Polished)


    def languageChange(self):
        self.setCaption(self.__tr("Paul Browser"))
        self.Dirlist.header().setLabel(0,self.__tr("Column 1"))
        self.Dirlist.clear()
        item = QListViewItem(self.Dirlist,None)
        item.setText(0,self.__tr("New Item"))



    def __tr(self,s,c = None):
        return qApp.translate("BrowserWindow",s,c)
