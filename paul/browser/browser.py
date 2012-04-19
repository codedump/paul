#!/usr/bin/python

from paul.browser.ui_browserwindow import *
from PyQt4 import QtGui
import sys

testfile = "/home/florin/local/analysis/uru2si2/2010-zpoint/jul2010.uxp-dir/jul10_urs11/t10k/jul10_urs11_09gif.ibw"

app = QtGui.QApplication (sys.argv)
win = QtGui.QMainWindow ()
win_ui = Ui_BrowserWindow()
win_ui.setupUi (win)

win.show()

sys.exit(app.exec_())
