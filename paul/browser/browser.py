#!/usr/bin/python

from paul.browser.browserwindow import BrowserWindow
from PyQt4 import QtGui
import sys

testfile = "/home/florin/local/analysis/uru2si2/2010-zpoint/jul2010.uxp-dir/jul10_urs11/t10k/jul10_urs11_09gif.ibw"

def run():
    app = QtGui.QApplication (sys.argv)
    win = BrowserWindow()
    win.show()

    sys.exit(app.exec_())

if __name__ == "__main__":
    run()
